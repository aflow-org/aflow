// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2019
#ifndef _AFLOW_XATOM_CPP
#define _AFLOW_XATOM_CPP

#include "structure/aflow_xatom.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <ios>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_automatic_template.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_defs.h"
#include "aflow_xhost.h"
#include "flow/aflow_ivasp.h"
#include "modules/SYM/aflow_symmetry_spacegroup.h" //DX20180723

// constructors
_atom::_atom() {
  free();
}

// destructor
_atom::~_atom() {
  free();
}

void _atom::free() {
  // PRIVATE //HE20220826 changed constructor to use free()
  fpos.clear();
  cpos.clear();
  corigin.clear();
  coord.clear();
  fpos_equation.clear(); // DX20180607 - symbolic math for atom positions
  cpos_equation.clear(); // DX20180607 - symbolic math for atom positions
  spin = 0.0;
  spin_is_given = false; // DX20170921 - magnetic sym
  noncoll_spin.clear(); // DX20171205 - magnetic sym (non-collinear)
  noncoll_spin_is_given = false; // DX20171205 - magnetic sym (non-collinear)
  mass = 0.0;
  type = 0;
  name = "";
  name_is_given = false;
  cleanname = "";
  info = 0; // (RHT)
  atomic_number = 0;
  //[CO20200130 - number->basis]number=0;
  sd = "";
  ijk.clear();
  isincell = false;
  basis = 0;
  reference = 0.0;
  ireference = 0;
  equivalent = -1;
  is_inequivalent = true;
  num_equivalents = 0;
  index_iatoms = 0;
  order_parameter_value = 0;
  order_parameter_atom = false;
  partial_occupation_value = 1.0;
  partial_occupation_flag = false;
  shell = 0;
  force.clear(); // CO20211107
  verbose = false;
  print_RHT = false; // CO20190405 //true; //CHANGE THIS BACK TO false WHEN DONE DEBUGGING  (RHT)
  print_cartesian = false;
}

void _atom::copy(const _atom& b) {
  // copy PRIVATE
  fpos = b.fpos;
  cpos = b.cpos;
  corigin = b.corigin;
  coord = b.coord;
  fpos_equation = b.fpos_equation; // DX20180607 - symbolic math for atom positions
  cpos_equation = b.cpos_equation; // DX20180607 - symbolic math for atom positions
  spin = b.spin;
  spin_is_given = b.spin_is_given; // DX20170921 - magnetic sym
  noncoll_spin = b.noncoll_spin; // DX20171205 - magnetic sym (non-collinear)
  noncoll_spin_is_given = b.noncoll_spin_is_given; // DX20171205 - magnetic sym (non-collinear)
  mass = b.mass;
  type = b.type;
  name = b.name;
  name_is_given = b.name_is_given;
  cleanname = b.cleanname;
  info = b.info; // (RHT)
  atomic_number = b.atomic_number;
  //[CO20200130 - number->basis]number=b.number;
  sd = b.sd;
  ijk = b.ijk;
  isincell = b.isincell;
  basis = b.basis;
  reference = b.reference;
  ireference = b.ireference;
  equivalent = b.equivalent;
  is_inequivalent = b.is_inequivalent;
  num_equivalents = b.num_equivalents;
  index_iatoms = b.index_iatoms;
  order_parameter_value = b.order_parameter_value;
  order_parameter_atom = b.order_parameter_atom;
  partial_occupation_value = b.partial_occupation_value;
  partial_occupation_flag = b.partial_occupation_flag;
  shell = b.shell;
  force = b.force; // CO20211107
  verbose = b.verbose;
  print_RHT = b.print_RHT; // (RHT)
  print_cartesian = b.print_cartesian;
}

const _atom& _atom::operator=(const _atom& b) {
  // operator= PUBLIC
  if (this != &b) {
    copy(b);
  }
  return *this;
}

_atom::_atom(const _atom& b) {
  // copy PUBLIC
  //  free(); *this=b;
  copy(b);
}

void _atom::clear() {
  free();
}

ostream& operator<<(ostream& oss, const _atom& atom) {
  oss.setf(std::ios::fixed, std::ios::floatfield);
  oss.precision(10);
  if (atom.print_RHT == true) {
    // oss << "ATOM COUT-RHT" << endl;
    oss << atom.coord << " " << atom.name;
  } else {
    if (atom.verbose == true) {
      oss << " " << endl;
      oss << "ATOM COUT" << endl;
      oss << "type=" << atom.type << endl;
      oss << "spin=" << atom.spin << endl;
      oss << "spin_is_given=" << atom.spin_is_given << endl; // DX20170921 - magnetic sym
      oss << "noncoll_spin=" << atom.noncoll_spin << endl; // DX20171205 - magnetic sym (non-collinear)
      oss << "noncoll_spin_is_given=" << atom.noncoll_spin_is_given << endl; // DX20171205 - magnetic sym (non-collinear)
      oss << "mass=" << atom.mass << endl;
      oss << "name=" << atom.name << endl;
      oss << "info=" << atom.info << endl;
      oss << "cleanname=" << atom.cleanname << endl;
      oss << "atomic_number=" << atom.atomic_number << endl;
      oss << "basis=" << atom.basis << endl; //[CO20200130 - number->basis]oss << "number=" << atom.number << endl;
      oss << "name_is_given=" << atom.name_is_given << endl;
      oss << "sd=" << atom.sd << endl;
      oss << "print_cartesian" << atom.print_cartesian << endl;
      oss << "print_RHT" << atom.print_RHT << endl;
      oss << "corigin" << atom.corigin(1) << " " << atom.corigin(2) << " " << atom.corigin(3) << endl;
      oss << "coord" << atom.coord(1) << " " << atom.coord(2) << " " << atom.coord(3) << endl;

      oss << "fpos_equation" << aurostd::joinWDelimiter(atom.fpos_equation, " ") << endl;
      // DX20180607 - symbolic math for atom positions //DX20191218 - join with delimiter in case empty
      oss << "cpos_equation" << aurostd::joinWDelimiter(atom.cpos_equation, " ") << endl;
      // DX20180607 - symbolic math for atom positions //DX20191218 - join with delimiter in case empty
      oss << "isincell=" << atom.isincell << endl;
      oss << "reference=" << atom.reference << endl;
      oss << "ireference=" << atom.ireference << endl;
      oss << "equivalent=" << atom.equivalent << endl;
      oss << "is_inequivalent=" << atom.is_inequivalent << endl;
      oss << "num_equivalents=" << atom.num_equivalents << endl;
      oss << "index_iatoms=" << atom.index_iatoms << endl;
      oss << "verbose=" << atom.verbose << endl;
      oss << "order_parameter_value=" << atom.order_parameter_value << endl;
      oss << "order_parameter_atom=" << atom.order_parameter_atom << endl;
      oss << "partial_occupation_value=" << atom.partial_occupation_value << endl;
      oss << "partial_occupation_flag=" << atom.partial_occupation_flag << endl;
      oss << "nearest_neighbor_shell_num= " << atom.shell << endl;
      oss << "force= " << atom.force << endl;
    }
    if (atom.print_cartesian == true) {
      if (atom.verbose) {
        oss << "cartesian" << endl;
      }
      oss << "C " << atom.cpos(1) << " " << atom.cpos(2) << " " << atom.cpos(3); // << endl;
    } else {
      if (atom.verbose) {
        oss << "fractional" << endl;
      }
      oss << "F " << atom.fpos(1) << " " << atom.fpos(2) << " " << atom.fpos(3); //  << endl;
    }
    oss << " T=" << atom.type;
    oss << " B=" << atom.basis;
    //[CO20200130 - number->basis]oss << " N=" << atom.number;
    //  oss << setw(1);
    oss << " ijk=[" << atom.ijk(1) << "," << atom.ijk(2) << "," << atom.ijk(3) << "]";
    if (atom.verbose) {
      oss << " " << endl;
    }
  }
  return oss;
}

void _atom::CleanName() {
  // CO20200624 - cleaning up function
  // the old function had a sign problem and was inefficient
  cleanname = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(name);
  aurostd::VASP_PseudoPotential_CleanName_InPlace(cleanname);
  // DX20200907 - changed from KBIN to aurostd (since this is now in xparser)
  if (cleanname.size() > 3) {
    cleanname = cleanname.substr(0, 3);
  } // cannot be longer than 3 characters
  for (size_t i = cleanname.size() - 1; i < cleanname.size(); i--) {
    // go backwards and clean anything that isn't between A-Z and a-z
    // 65-90 is A-Z
    // 97-122 is a-z
    if (!((cleanname[i] >= 65 && cleanname[i] <= 90) || (cleanname[i] >= 97 && cleanname[i] <= 122))) {
      cleanname.erase(cleanname.begin() + i);
    }
  }
  // the old function insists on correcting if we have aL instead of Al...
  if (!cleanname.empty() && cleanname[0] >= 97 && cleanname[0] <= 122) {
    cleanname[0] -= ('a' - 'A');
  }
  for (size_t i = 1; i < cleanname.size(); i++) {
    if (cleanname[i] >= 65 && cleanname[i] <= 90) {
      cleanname[i] += ('a' - 'A');
    }
  }

  for (size_t j = 0; j < vatom_symbol.size(); j++) {
    if (cleanname == vatom_symbol[j]) {
      atomic_number = j;
    }
  }
}

void _atom::CleanSpin() {
  spin = 0.0;
  spin_is_given = false; // DX20170921 - magnetic sym
  noncoll_spin.clear(); // DX20171205 - magnetic sym (non-collinear)
  noncoll_spin_is_given = false; // DX20171205 - magnetic sym (non-collinear)
  if (!XHOST.READ_SPIN_FROM_ATOMLABEL) {
    return;
  } // SD20220316
  if (name.find("+") != string::npos) {
    spin = atof(name.substr(name.find("+")).c_str());
    spin_is_given = true;
  } // DX20170921 - magnetic sym
  if (name.find("-") != string::npos) {
    spin = atof(name.substr(name.find("-")).c_str());
    spin_is_given = true;
  } // DX20170921 - magnetic sym
}

void _atom::ClearSymmetry() {
  // CO20190219
  (*this).equivalent = -1;
  (*this).is_inequivalent = true;
  (*this).num_equivalents = 0;
}

std::vector<string> vatom_symbol(NUM_ELEMENTS); // store starting from ONE
std::vector<string> vatom_name(NUM_ELEMENTS); // store starting from ONE
std::vector<double> vatom_mass(NUM_ELEMENTS); // store starting from ONE
std::vector<double> vatom_volume(NUM_ELEMENTS); // store starting from ONE
std::vector<int> vatom_valence_iupac(NUM_ELEMENTS);
// store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
std::vector<int> vatom_valence_std(NUM_ELEMENTS);
// store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
std::vector<double> vatom_miedema_phi_star(NUM_ELEMENTS);
// store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> vatom_miedema_nws(NUM_ELEMENTS);
// store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> vatom_miedema_Vm(NUM_ELEMENTS);
// store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> vatom_miedema_gamma_s(NUM_ELEMENTS);
// store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
std::vector<double> vatom_miedema_BVm(NUM_ELEMENTS);
// store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
// for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
std::vector<double> vatom_radius(NUM_ELEMENTS); // store starting from ONE
std::vector<double> vatom_radius_covalent(NUM_ELEMENTS); // store starting from ONE //DX+CO20170904
std::vector<double> vatom_electronegativity(NUM_ELEMENTS); // store starting from ONE
std::vector<string> vatom_crystal(NUM_ELEMENTS); // store starting from ONE
std::vector<double> vatom_xray_scatt(NUM_ELEMENTS); // store starting from ONE
std::vector<double> vatom_pettifor_scale(NUM_ELEMENTS);
// store starting from ONE Chemical Scale Pettifor Solid State Communications 51 31-34 1984
std::vector<double> vatom_pearson_coefficient(NUM_ELEMENTS); // ME20181020 Pearson mass deviation coefficient

void atoms_initialize() {
  for (int i = 0; i < NUM_ELEMENTS; i++) {
    // clear
    vatom_xray_scatt.at(i) = (double) i - 1; // shift+1
    vatom_mass.at(i) = AMU2KILOGRAM; // masses in kilos
    vatom_volume.at(i) = NNN; // atomic volume in A^3 from the FCC vasp table and/or successive calculations
    vatom_valence_iupac.at(i) = NNN;
    // IUPAC Maximum number of univalent atoms that may combine with an atom of the element under consideration, or with a fragment, or for which an atom of this element can be substituted.
    vatom_valence_std.at(i) = NNN; // stanmdard: number electrons minus closed shell at leff (noble gas)
    vatom_miedema_phi_star.at(i) = NNN; // Miedema Rule Table 1a Physica 100B (1980) 1-28   (phi^\star in (V))
    vatom_miedema_nws.at(i) = NNN; // Miedema Rule Table 1a Physica 100B (1980) 1-28   n_{ws}^{1/3} in (d.u.)^1/3
    vatom_miedema_Vm.at(i) = NNN; // Miedema Rule Table 1a Physica 100B (1980) 1-28   V_m^{2/3} in (cm^2)
    vatom_miedema_gamma_s.at(i) = NNN; // Miedema Rule Table 1a Physica 100B (1980) 1-28   \gamma_s^0 in (mJ/m^2)
    vatom_miedema_BVm.at(i) = NNN; // Miedema Rule Table 1a Physica 100B (1980) 1-28   BV_m (kJ/mole)
    vatom_radius.at(i) = NNN; // Saxena (nm)
    vatom_radius_covalent.at(i) = NNN; // Codero (Angstroms) //DX+CO20170904
    vatom_electronegativity.at(i) = NNN; // Saxena
    vatom_crystal.at(i) = "nnn"; // Ashcroft-Mermin
  }

  // Xray_scatt_vector All data collected from the NIST online tables
  // http://physics.nist.gov/PhysRefData/FFast/html/form.html
  // All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  // These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).

  // All data collected from the online tables:
  // http://www-cxro.lbl.gov/optical_constants/pert_form.html
  // All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV].

  int i;
  // ROW 0
  i = 0;
  vatom_symbol[i] = "XX";
  vatom_name[i] = "UNDEFINED";
  vatom_mass[i] *= 0.00;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 0;
  vatom_valence_iupac[i] = 0;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = NNN;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_xray_scatt[i] = 0;
  vatom_pearson_coefficient[i] = 0.0;

  // ROW 1
  // s-electron systems
  i++;
  vatom_symbol[i] = "H";
  vatom_name[i] = "Hydrogen";
  vatom_mass[i] *= 1.0079;
  vatom_volume[i] = 0.75110;
  vatom_valence_std[i] = 1;
  vatom_valence_iupac[i] = 1;
  vatom_miedema_phi_star[i] = 5.2;
  vatom_miedema_nws[i] = 1.5;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.046;
  vatom_radius_covalent[i] = 0.31;
  vatom_electronegativity[i] = 2.10;
  vatom_crystal[i] = "hex";
  vatom_pettifor_scale[i] = 0;
  vatom_xray_scatt[i] = 1.000;
  vatom_pearson_coefficient[i] = 0.00011460743; // H volume wrong *dimer*   MIEDEMA = PAUL VAN DER PUT book
  i++;
  vatom_symbol[i] = "He";
  vatom_name[i] = "Helium";
  vatom_mass[i] *= 4.0026;
  vatom_volume[i] = -1.000;
  vatom_valence_std[i] = 0;
  vatom_valence_iupac[i] = 0;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 0.28;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0;
  vatom_xray_scatt[i] = 2.000;
  vatom_pearson_coefficient[i] = 8.32328E-8; // He

  // ROW2
  // s-electron systems
  i++;
  vatom_symbol[i] = "Li";
  vatom_name[i] = "Lithium";
  vatom_mass[i] *= 6.941;
  vatom_volume[i] = 20.24110;
  vatom_valence_std[i] = 1;
  vatom_valence_iupac[i] = 1;
  vatom_miedema_phi_star[i] = 2.85;
  vatom_miedema_nws[i] = 0.98;
  vatom_miedema_Vm[i] = 5.5;
  vatom_miedema_gamma_s[i] = 530;
  vatom_miedema_BVm[i] = 1.5;
  vatom_radius[i] = 0.152;
  vatom_radius_covalent[i] = 1.28;
  vatom_electronegativity[i] = 0.98;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.45;
  vatom_xray_scatt[i] = 3.00145;
  vatom_pearson_coefficient[i] = 0.0014588232; // Li
  i++;
  vatom_symbol[i] = "Be";
  vatom_name[i] = "Beryllium";
  vatom_mass[i] *= 9.0122;
  vatom_volume[i] = 7.83290;
  vatom_valence_std[i] = 2;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = 4.20;
  vatom_miedema_nws[i] = 1.60;
  vatom_miedema_Vm[i] = 2.9;
  vatom_miedema_gamma_s[i] = 1900;
  vatom_miedema_BVm[i] = 4.9;
  vatom_radius[i] = 0.114;
  vatom_radius_covalent[i] = 0.96;
  vatom_electronegativity[i] = 1.57;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 1.50;
  vatom_pearson_coefficient[i] = 0.0; // Be
  // p-electron systems
  i++;
  vatom_symbol[i] = "B";
  vatom_name[i] = "Boron";
  vatom_mass[i] *= 10.81;
  vatom_volume[i] = 5.88420;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 4.75;
  vatom_miedema_nws[i] = 1.55;
  vatom_miedema_Vm[i] = 2.8;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.097;
  vatom_radius_covalent[i] = 0.84;
  vatom_electronegativity[i] = 2.04;
  vatom_crystal[i] = "tet";
  vatom_pettifor_scale[i] = 2.00;
  vatom_pearson_coefficient[i] = 0.00135391428; // B
  i++;
  vatom_symbol[i] = "C";
  vatom_name[i] = "Carbon";
  vatom_mass[i] *= 12.011;
  vatom_volume[i] = 5.59490;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 6.20;
  vatom_miedema_nws[i] = 1.90;
  vatom_miedema_Vm[i] = 1.8;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.077;
  vatom_radius_covalent[i] = 0.76;
  vatom_electronegativity[i] = 2.55;
  vatom_crystal[i] = "dia";
  vatom_pettifor_scale[i] = 2.50;
  vatom_xray_scatt[i] = 6.019;
  vatom_pearson_coefficient[i] = 0.00007387218;
  // C  //DX+CO20170904 vatom_radius_covalent uses sp3 hybridization (most common)
  i++;
  vatom_symbol[i] = "N";
  vatom_name[i] = "Nitrogen";
  vatom_mass[i] *= 14.0067;
  vatom_volume[i] = 7.59940;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 7.00;
  vatom_miedema_nws[i] = 1.60;
  vatom_miedema_Vm[i] = 2.2;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.071;
  vatom_radius_covalent[i] = 0.71;
  vatom_electronegativity[i] = 3.04;
  vatom_crystal[i] = "hex";
  vatom_pettifor_scale[i] = 3.00;
  vatom_pearson_coefficient[i] = 0.00001857771; // N JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "O";
  vatom_name[i] = "Oxygen";
  vatom_mass[i] *= 15.9994;
  vatom_volume[i] = 7.78230;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = 6.97;
  vatom_miedema_nws[i] = 1.70;
  vatom_miedema_Vm[i] = 2.656;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.060;
  vatom_radius_covalent[i] = 0.66;
  vatom_electronegativity[i] = 3.44;
  vatom_crystal[i] = "cub";
  vatom_pettifor_scale[i] = 3.50;
  vatom_xray_scatt[i] = 8.052;
  vatom_pearson_coefficient[i] = 0.00003358805; // O Table 27 of JX
  i++;
  vatom_symbol[i] = "F";
  vatom_name[i] = "Fluorine";
  vatom_mass[i] *= 18.9984;
  vatom_volume[i] = 9.99090;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 1;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 0.57;
  vatom_electronegativity[i] = 3.98;
  vatom_crystal[i] = "mcl";
  vatom_pettifor_scale[i] = 4.00;
  vatom_pearson_coefficient[i] = 0.0; // F
  i++;
  vatom_symbol[i] = "Ne";
  vatom_name[i] = "Neon";
  vatom_mass[i] *= 20.179;
  vatom_volume[i] = 19.9052;
  vatom_valence_std[i] = 0;
  vatom_valence_iupac[i] = 0;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.160;
  vatom_radius_covalent[i] = 0.58;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.00082783369; // Ne volume calculated with fcc-pawpbe

  // ROW3
  // s-electron systems
  i++;
  vatom_symbol[i] = "Na";
  vatom_name[i] = "Sodium";
  vatom_mass[i] *= 22.9898;
  vatom_volume[i] = 36.9135;
  vatom_valence_std[i] = 1;
  vatom_valence_iupac[i] = 1;
  vatom_miedema_phi_star[i] = 2.70;
  vatom_miedema_nws[i] = 0.82;
  vatom_miedema_Vm[i] = 8.3;
  vatom_miedema_gamma_s[i] = 260;
  vatom_miedema_BVm[i] = 1.6;
  vatom_radius[i] = 0.186;
  vatom_radius_covalent[i] = 1.66;
  vatom_electronegativity[i] = 0.93;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.40;
  vatom_pearson_coefficient[i] = 0.0; // Na
  i++;
  vatom_symbol[i] = "Mg";
  vatom_name[i] = "Magnesium";
  vatom_mass[i] *= 24.305;
  vatom_volume[i] = 22.8178;
  vatom_valence_std[i] = 2;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = 3.45;
  vatom_miedema_nws[i] = 1.17;
  vatom_miedema_Vm[i] = 5.8;
  vatom_miedema_gamma_s[i] = 790;
  vatom_miedema_BVm[i] = 5.0;
  vatom_radius[i] = 0.160;
  vatom_radius_covalent[i] = 1.41;
  vatom_electronegativity[i] = 1.31;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 1.28;
  vatom_pearson_coefficient[i] = 0.00073988271; // Mg
  // p-electron systems
  i++;
  vatom_symbol[i] = "Al";
  vatom_name[i] = "Aluminium";
  vatom_mass[i] *= 26.9815;
  vatom_volume[i] = 16.4000;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 4.20;
  vatom_miedema_nws[i] = 1.39;
  vatom_miedema_Vm[i] = 4.6;
  vatom_miedema_gamma_s[i] = 1200;
  vatom_miedema_BVm[i] = 7.2;
  vatom_radius[i] = 0.143;
  vatom_radius_covalent[i] = 1.21;
  vatom_electronegativity[i] = 1.61;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.66;
  vatom_pearson_coefficient[i] = 0.0; // Al
  i++;
  vatom_symbol[i] = "Si";
  vatom_name[i] = "Silicon";
  vatom_mass[i] *= 28.0855;
  vatom_volume[i] = 14.3536;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 4.70;
  vatom_miedema_nws[i] = 1.50;
  vatom_miedema_Vm[i] = 4.2;
  vatom_miedema_gamma_s[i] = 1290;
  vatom_miedema_BVm[i] = 11.9;
  vatom_radius[i] = 0.117;
  vatom_radius_covalent[i] = 1.11;
  vatom_electronegativity[i] = 1.90;
  vatom_crystal[i] = "dia";
  vatom_pettifor_scale[i] = 1.92;
  vatom_xray_scatt[i] = 14.43;
  vatom_pearson_coefficient[i] = 0.00020046752; // Si ???
  i++;
  vatom_symbol[i] = "P";
  vatom_name[i] = "Phosphorus";
  vatom_mass[i] *= 30.9738;
  vatom_volume[i] = 14.1995;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 5.5;
  vatom_miedema_nws[i] = 1.65;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.109;
  vatom_radius_covalent[i] = 1.07;
  vatom_electronegativity[i] = 2.19;
  vatom_crystal[i] = "cub";
  vatom_pettifor_scale[i] = 2.18;
  vatom_xray_scatt[i] = 15.3133;
  vatom_pearson_coefficient[i] = 0.0; // P   MIEDEMA = PAUL VAN DER PUT book
  i++;
  vatom_symbol[i] = "S";
  vatom_name[i] = "Sulphur";
  vatom_mass[i] *= 32.06;
  vatom_volume[i] = 15.7301;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 5.6;
  vatom_miedema_nws[i] = 1.46;
  vatom_miedema_Vm[i] = 4.376;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.106;
  vatom_radius_covalent[i] = 1.05;
  vatom_electronegativity[i] = 2.58;
  vatom_crystal[i] = "orc";
  vatom_pettifor_scale[i] = 2.44;
  vatom_pearson_coefficient[i] = 0.00016807795; // S Table 27 of JX
  i++;
  vatom_symbol[i] = "Cl";
  vatom_name[i] = "Chlorine";
  vatom_mass[i] *= 35.453;
  vatom_volume[i] = 21.2947;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = 5.32;
  vatom_miedema_nws[i] = 0.34;
  vatom_miedema_Vm[i] = 6.71;
  vatom_miedema_gamma_s[i] = 1013;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.107;
  vatom_radius_covalent[i] = 1.02;
  vatom_electronegativity[i] = 3.16;
  vatom_crystal[i] = "orc";
  vatom_pettifor_scale[i] = 2.70;
  vatom_pearson_coefficient[i] = 0.00058238731; // Cl interpolation phi_star, nws, Vm, gamma  JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Ar";
  vatom_name[i] = "Argon";
  vatom_mass[i] *= 39.948;
  vatom_volume[i] = 22.000;
  vatom_valence_std[i] = 0;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.192;
  vatom_radius_covalent[i] = 1.06;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.00003509919;
  // Ar guessed volume, must double check from results  JX CHANGED VALENCE

  // ROW4
  // s-electron systems
  i++;
  vatom_symbol[i] = "K";
  vatom_name[i] = "Potassium";
  vatom_mass[i] *= 39.0983;
  vatom_volume[i] = 73.9091;
  vatom_valence_std[i] = 1;
  vatom_valence_iupac[i] = 1;
  vatom_miedema_phi_star[i] = 2.25;
  vatom_miedema_nws[i] = 0.65;
  vatom_miedema_Vm[i] = 12.8;
  vatom_miedema_gamma_s[i] = 150;
  vatom_miedema_BVm[i] = 1.5;
  vatom_radius[i] = 0.231;
  vatom_radius_covalent[i] = 2.03;
  vatom_electronegativity[i] = 0.82;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0.35;
  vatom_pearson_coefficient[i] = 0.000164; // K
  i++;
  vatom_symbol[i] = "Ca";
  vatom_name[i] = "Calcium";
  vatom_mass[i] *= 40.08;
  vatom_volume[i] = 42.1927;
  vatom_valence_std[i] = 2;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = 2.55;
  vatom_miedema_nws[i] = 0.91;
  vatom_miedema_Vm[i] = 8.8;
  vatom_miedema_gamma_s[i] = 490;
  vatom_miedema_BVm[i] = 4.0;
  vatom_radius[i] = 0.197;
  vatom_radius_covalent[i] = 1.76;
  vatom_electronegativity[i] = 1.00;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.60;
  vatom_pearson_coefficient[i] = 0.000297564; // Ca
  // d-electron systems: transition metals
  i++;
  vatom_symbol[i] = "Sc";
  vatom_name[i] = "Scandium";
  vatom_mass[i] *= 44.9559;
  vatom_volume[i] = 24.6739;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.25;
  vatom_miedema_nws[i] = 1.27;
  vatom_miedema_Vm[i] = 6.1;
  vatom_miedema_gamma_s[i] = 1200;
  vatom_miedema_BVm[i] = 6.6;
  vatom_radius[i] = 0.160;
  vatom_radius_covalent[i] = 1.70;
  vatom_electronegativity[i] = 1.36;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.74;
  vatom_xray_scatt[i] = 21.34;
  vatom_pearson_coefficient[i] = 0.0; // Sc
  i++;
  vatom_symbol[i] = "Ti";
  vatom_name[i] = "Titanium";
  vatom_mass[i] *= 47.9;
  vatom_volume[i] = 17.1035;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.65;
  vatom_miedema_nws[i] = 1.47;
  vatom_miedema_Vm[i] = 4.8;
  vatom_miedema_gamma_s[i] = 2050;
  vatom_miedema_BVm[i] = 11.0;
  vatom_radius[i] = 0.147;
  vatom_radius_covalent[i] = 1.60;
  vatom_electronegativity[i] = 1.54;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.79;
  vatom_xray_scatt[i] = 22.24;
  vatom_pearson_coefficient[i] = 0.000286456; // Ti
  i++;
  vatom_symbol[i] = "V";
  vatom_name[i] = "Vanadium";
  vatom_mass[i] *= 50.9415;
  vatom_volume[i] = 13.2086;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 4.25;
  vatom_miedema_nws[i] = 1.64;
  vatom_miedema_Vm[i] = 4.1;
  vatom_miedema_gamma_s[i] = 2600;
  vatom_miedema_BVm[i] = 14.0;
  vatom_radius[i] = 0.132;
  vatom_radius_covalent[i] = 1.53;
  vatom_electronegativity[i] = 1.63;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.84;
  vatom_pearson_coefficient[i] = 9.54831E-07; // V
  i++;
  vatom_symbol[i] = "Cr";
  vatom_name[i] = "Chromium";
  vatom_mass[i] *= 51.996;
  vatom_volume[i] = 11.4136;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 4.65;
  vatom_miedema_nws[i] = 1.74;
  vatom_miedema_Vm[i] = 3.7;
  vatom_miedema_gamma_s[i] = 2400;
  vatom_miedema_BVm[i] = 14.0;
  vatom_radius[i] = 0.125;
  vatom_radius_covalent[i] = 1.39;
  vatom_electronegativity[i] = 1.66;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.89;
  vatom_xray_scatt[i] = 23.84;
  vatom_pearson_coefficient[i] = 0.00013287; // Cr
  i++;
  vatom_symbol[i] = "Mn";
  vatom_name[i] = "Manganese";
  vatom_mass[i] *= 54.93805;
  vatom_volume[i] = 10.6487;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = 4.45;
  vatom_miedema_nws[i] = 1.61;
  vatom_miedema_Vm[i] = 3.8;
  vatom_miedema_gamma_s[i] = 1600;
  vatom_miedema_BVm[i] = 4.4;
  vatom_radius[i] = 0.112;
  vatom_radius_covalent[i] = 1.61;
  vatom_electronegativity[i] = 1.55;
  vatom_crystal[i] = "cub";
  vatom_pettifor_scale[i] = 0.94;
  vatom_xray_scatt[i] = 24.46;
  vatom_pearson_coefficient[i] = 1.67276E-32;
  // vatom_xray_scatt[i]=24.3589; Mn JX CHANGED VALENCE //DX+CO20170904 vatom_radius_covalent[i] uses high spin configuration (most frequent)
  i++;
  vatom_symbol[i] = "Fe";
  vatom_name[i] = "Iron";
  vatom_mass[i] *= 55.847;
  vatom_volume[i] = 10.2315;
  vatom_valence_std[i] = 8;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 4.93;
  vatom_miedema_nws[i] = 1.77;
  vatom_miedema_Vm[i] = 3.7;
  vatom_miedema_gamma_s[i] = 2550;
  vatom_miedema_BVm[i] = 12.0;
  vatom_radius[i] = 0.124;
  vatom_radius_covalent[i] = 1.52;
  vatom_electronegativity[i] = 1.83;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.99;
  vatom_xray_scatt[i] = 24.85;
  vatom_pearson_coefficient[i] = 9.17912E-05;
  // vatom_xray_scatt[i]=24.6830; Fe JX CHANGED VALENCE //DX+CO20170904 vatom_radius_covalent[i] uses high spin configuration (most frequent)
  i++;
  vatom_symbol[i] = "Co";
  vatom_name[i] = "Cobalt";
  vatom_mass[i] *= 58.9332;
  vatom_volume[i] = 10.3205;
  vatom_valence_std[i] = 9;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 5.10;
  vatom_miedema_nws[i] = 1.75;
  vatom_miedema_Vm[i] = 3.5;
  vatom_miedema_gamma_s[i] = 2550;
  vatom_miedema_BVm[i] = 13.0;
  vatom_radius[i] = 0.125;
  vatom_radius_covalent[i] = 1.26;
  vatom_electronegativity[i] = 1.88;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 1.04;
  vatom_xray_scatt[i] = 24.59;
  vatom_pearson_coefficient[i] = 0.0;
  // Co JX CHANGED VALENCE //DX+CO20170904 vatom_radius_covalent[i] uses low spin configuration (most frequent)
  i++;
  vatom_symbol[i] = "Ni";
  vatom_name[i] = "Nickel";
  vatom_mass[i] *= 58.69;
  vatom_volume[i] = 10.8664;
  vatom_valence_std[i] = 10;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 5.20;
  vatom_miedema_nws[i] = 1.75;
  vatom_miedema_Vm[i] = 3.5;
  vatom_miedema_gamma_s[i] = 2450;
  vatom_miedema_BVm[i] = 12.0;
  vatom_radius[i] = 0.125;
  vatom_radius_covalent[i] = 1.24;
  vatom_electronegativity[i] = 1.91;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.09;
  vatom_xray_scatt[i] = 25.02;
  vatom_pearson_coefficient[i] = 0.000430773; // Ni
  i++;
  vatom_symbol[i] = "Cu";
  vatom_name[i] = "Copper";
  vatom_mass[i] *= 63.546;
  vatom_volume[i] = 12.0159;
  vatom_valence_std[i] = 11;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 4.55;
  vatom_miedema_nws[i] = 1.47;
  vatom_miedema_Vm[i] = 3.7;
  vatom_miedema_gamma_s[i] = 1850;
  vatom_miedema_BVm[i] = 9.3;
  vatom_radius[i] = 0.128;
  vatom_radius_covalent[i] = 1.32;
  vatom_electronegativity[i] = 1.90;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.20;
  vatom_xray_scatt[i] = 27.03;
  vatom_pearson_coefficient[i] = 0.00021086; // Cu JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Zn";
  vatom_name[i] = "Zinc";
  vatom_mass[i] *= 65.38;
  vatom_volume[i] = 15.0827;
  vatom_valence_std[i] = 12;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = 4.10;
  vatom_miedema_nws[i] = 1.32;
  vatom_miedema_Vm[i] = 4.4;
  vatom_miedema_gamma_s[i] = 1020;
  vatom_miedema_BVm[i] = 5.5;
  vatom_radius[i] = 0.133;
  vatom_radius_covalent[i] = 1.22;
  vatom_electronegativity[i] = 1.65;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 1.44;
  vatom_xray_scatt[i] = 28.44;
  vatom_pearson_coefficient[i] = 0.000595597; // Zn
  // p-electron systems
  i++;
  vatom_symbol[i] = "Ga";
  vatom_name[i] = "Gallium";
  vatom_mass[i] *= 69.737;
  vatom_volume[i] = 18.9039;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 4.10;
  vatom_miedema_nws[i] = 1.31;
  vatom_miedema_Vm[i] = 5.2;
  vatom_miedema_gamma_s[i] = 830;
  vatom_miedema_BVm[i] = 6.7;
  vatom_radius[i] = 0.135;
  vatom_radius_covalent[i] = 1.22;
  vatom_electronegativity[i] = 1.81;
  vatom_crystal[i] = "orc";
  vatom_pettifor_scale[i] = 1.68;
  vatom_pearson_coefficient[i] = 0.000197588; // Ga
  i++;
  vatom_symbol[i] = "Ge";
  vatom_name[i] = "Germanium";
  vatom_mass[i] *= 72.59;
  vatom_volume[i] = 19.2948;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 4.55;
  vatom_miedema_nws[i] = 1.37;
  vatom_miedema_Vm[i] = 4.6;
  vatom_miedema_gamma_s[i] = 1030;
  vatom_miedema_BVm[i] = 10.5;
  vatom_radius[i] = 0.122;
  vatom_radius_covalent[i] = 1.20;
  vatom_electronegativity[i] = 2.01;
  vatom_crystal[i] = "dia";
  vatom_pettifor_scale[i] = 1.92;
  vatom_pearson_coefficient[i] = 0.00058782; // Ge
  i++;
  vatom_symbol[i] = "As";
  vatom_name[i] = "Arsenic";
  vatom_mass[i] *= 74.9216;
  vatom_volume[i] = 19.0677;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 4.80;
  vatom_miedema_nws[i] = 1.44;
  vatom_miedema_Vm[i] = 5.2;
  vatom_miedema_gamma_s[i] = 1000;
  vatom_miedema_BVm[i] = 5.1;
  vatom_radius[i] = 0.125;
  vatom_radius_covalent[i] = 1.19;
  vatom_electronegativity[i] = 2.18;
  vatom_crystal[i] = "rhl";
  vatom_pettifor_scale[i] = 2.16;
  vatom_pearson_coefficient[i] = 0.0; // As
  i++;
  vatom_symbol[i] = "Se";
  vatom_name[i] = "Selenium";
  vatom_mass[i] *= 78.96;
  vatom_volume[i] = 20.3733;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 5.17;
  vatom_miedema_nws[i] = 1.40;
  vatom_miedema_Vm[i] = 5.172;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.116;
  vatom_radius_covalent[i] = 1.20;
  vatom_electronegativity[i] = 2.55;
  vatom_crystal[i] = "hex";
  vatom_pettifor_scale[i] = 2.40;
  vatom_pearson_coefficient[i] = 0.00046279; // Se Table 27 of JX
  i++;
  vatom_symbol[i] = "Br";
  vatom_name[i] = "Bromine";
  vatom_mass[i] *= 79.904;
  vatom_volume[i] = 26.3292;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = 5.20;
  vatom_miedema_nws[i] = 1.35;
  vatom_miedema_Vm[i] = 7.31;
  vatom_miedema_gamma_s[i] = 943;
  vatom_miedema_BVm[i] = 3.4;
  vatom_radius[i] = 0.119;
  vatom_radius_covalent[i] = 1.20;
  vatom_electronegativity[i] = 2.96;
  vatom_crystal[i] = "orc";
  vatom_pettifor_scale[i] = 2.64;
  vatom_pearson_coefficient[i] = 0.000156277; // Br interpolation phi_star, nws, Vm, gamma, BVm JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Kr";
  vatom_name[i] = "Krypton";
  vatom_mass[i] *= 83.8;
  vatom_volume[i] = -1.0000;
  vatom_valence_std[i] = 0;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.197;
  vatom_radius_covalent[i] = 1.16;
  vatom_electronegativity[i] = 3;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.000248482; // Kr

  // ROW5
  // s-electron systems
  i++;
  vatom_symbol[i] = "Rb";
  vatom_name[i] = "Rubidium";
  vatom_mass[i] *= 85.4678;
  vatom_volume[i] = 91.2738;
  vatom_valence_std[i] = 1;
  vatom_valence_iupac[i] = 1;
  vatom_miedema_phi_star[i] = 2.10;
  vatom_miedema_nws[i] = 0.60;
  vatom_miedema_Vm[i] = 14.6;
  vatom_miedema_gamma_s[i] = 120;
  vatom_miedema_BVm[i] = 1.8;
  vatom_radius[i] = 0.251;
  vatom_radius_covalent[i] = 2.20;
  vatom_electronegativity[i] = 0.82;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.30;
  vatom_pearson_coefficient[i] = 0.000109697; // Rb
  i++;
  vatom_symbol[i] = "Sr";
  vatom_name[i] = "Strontium";
  vatom_mass[i] *= 87.62;
  vatom_volume[i] = 55.4105;
  vatom_valence_std[i] = 2;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = 2.40;
  vatom_miedema_nws[i] = 0.84;
  vatom_miedema_Vm[i] = 10.2;
  vatom_miedema_gamma_s[i] = 430;
  vatom_miedema_BVm[i] = 3.9;
  vatom_radius[i] = 0.215;
  vatom_radius_covalent[i] = 1.95;
  vatom_electronegativity[i] = 0.95;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0.55;
  vatom_pearson_coefficient[i] = 6.09969E-05; // Sr
  // d-electron systems: transition metals
  i++;
  vatom_symbol[i] = "Y";
  vatom_name[i] = "Yttrium";
  vatom_mass[i] *= 88.9059;
  vatom_volume[i] = 32.4546;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.20;
  vatom_miedema_nws[i] = 1.21;
  vatom_miedema_Vm[i] = 7.3;
  vatom_miedema_gamma_s[i] = 1100;
  vatom_miedema_BVm[i] = 7.2;
  vatom_radius[i] = 0.181;
  vatom_radius_covalent[i] = 1.90;
  vatom_electronegativity[i] = 1.22;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.70;
  vatom_pearson_coefficient[i] = 0.0; // Y
  i++;
  vatom_symbol[i] = "Zr";
  vatom_name[i] = "Zirconium";
  vatom_mass[i] *= 91.22;
  vatom_volume[i] = 23.2561;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.40;
  vatom_miedema_nws[i] = 1.39;
  vatom_miedema_Vm[i] = 5.8;
  vatom_miedema_gamma_s[i] = 1950;
  vatom_miedema_BVm[i] = 12.0;
  vatom_radius[i] = 0.158;
  vatom_radius_covalent[i] = 1.75;
  vatom_electronegativity[i] = 1.33;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.76;
  vatom_pearson_coefficient[i] = 0.000342629; // Zr
  i++;
  vatom_symbol[i] = "Nb";
  vatom_name[i] = "Niobium";
  vatom_mass[i] *= 92.9064;
  vatom_volume[i] = 18.3132;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 4.00;
  vatom_miedema_nws[i] = 1.62;
  vatom_miedema_Vm[i] = 4.9;
  vatom_miedema_gamma_s[i] = 2700;
  vatom_miedema_BVm[i] = 18.0;
  vatom_radius[i] = 0.143;
  vatom_radius_covalent[i] = 1.64;
  vatom_electronegativity[i] = 1.60;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.82;
  vatom_pearson_coefficient[i] = 0.0; // Nb
  i++;
  vatom_symbol[i] = "Mo";
  vatom_name[i] = "Molybdenum";
  vatom_mass[i] *= 95.94;
  vatom_volume[i] = 15.6175;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 4.65;
  vatom_miedema_nws[i] = 1.77;
  vatom_miedema_Vm[i] = 4.4;
  vatom_miedema_gamma_s[i] = 2950;
  vatom_miedema_BVm[i] = 26.0;
  vatom_radius[i] = 0.136;
  vatom_radius_covalent[i] = 1.54;
  vatom_electronegativity[i] = 2.16;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.88;
  vatom_pearson_coefficient[i] = 0.000598128; // Mo
  i++;
  vatom_symbol[i] = "Tc";
  vatom_name[i] = "Technetium";
  vatom_mass[i] *= 98.9062;
  vatom_volume[i] = 14.4670;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = 5.30;
  vatom_miedema_nws[i] = 1.81;
  vatom_miedema_Vm[i] = 4.2;
  vatom_miedema_gamma_s[i] = 3050;
  vatom_miedema_BVm[i] = 26.0;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 1.47;
  vatom_electronegativity[i] = 1.90;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.94;
  vatom_pearson_coefficient[i] = 0.0; // Tc JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Ru";
  vatom_name[i] = "Ruthenium";
  vatom_mass[i] *= 101.07;
  vatom_volume[i] = 13.8390;
  vatom_valence_std[i] = 8;
  vatom_valence_iupac[i] = 8;
  vatom_miedema_phi_star[i] = 5.40;
  vatom_miedema_nws[i] = 1.83;
  vatom_miedema_Vm[i] = 4.1;
  vatom_miedema_gamma_s[i] = 3050;
  vatom_miedema_BVm[i] = 26.0;
  vatom_radius[i] = 0.134;
  vatom_radius_covalent[i] = 1.46;
  vatom_electronegativity[i] = 2.20;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 1.00;
  vatom_pearson_coefficient[i] = 0.000406665; // Ru JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Rh";
  vatom_name[i] = "Rhodium";
  vatom_mass[i] *= 102.9055;
  vatom_volume[i] = 14.1731;
  vatom_valence_std[i] = 9;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 5.40;
  vatom_miedema_nws[i] = 1.76;
  vatom_miedema_Vm[i] = 4.1;
  vatom_miedema_gamma_s[i] = 2750;
  vatom_miedema_BVm[i] = 23.0;
  vatom_radius[i] = 0.134;
  vatom_radius_covalent[i] = 1.42;
  vatom_electronegativity[i] = 2.28;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.06;
  vatom_pearson_coefficient[i] = 1.90706E-32; // Rh
  i++;
  vatom_symbol[i] = "Pd";
  vatom_name[i] = "Palladium";
  vatom_mass[i] *= 106.4;
  vatom_volume[i] = 15.4596;
  vatom_valence_std[i] = 10;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 5.45;
  vatom_miedema_nws[i] = 1.67;
  vatom_miedema_Vm[i] = 4.3;
  vatom_miedema_gamma_s[i] = 2100;
  vatom_miedema_BVm[i] = 16.0;
  vatom_radius[i] = 0.137;
  vatom_radius_covalent[i] = 1.39;
  vatom_electronegativity[i] = 2.20;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.12;
  vatom_pearson_coefficient[i] = 0.000309478; // Pd
  i++;
  vatom_symbol[i] = "Ag";
  vatom_name[i] = "Silver";
  vatom_mass[i] *= 107.8682;
  vatom_volume[i] = 18.0678;
  vatom_valence_std[i] = 11;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 4.45;
  vatom_miedema_nws[i] = 1.39;
  vatom_miedema_Vm[i] = 4.7;
  vatom_miedema_gamma_s[i] = 1250;
  vatom_miedema_BVm[i] = 10.0;
  vatom_radius[i] = 0.144;
  vatom_radius_covalent[i] = 1.45;
  vatom_electronegativity[i] = 1.93;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.18;
  vatom_xray_scatt[i] = 47.18;
  vatom_pearson_coefficient[i] = 8.57985E-05; // Ag JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Cd";
  vatom_name[i] = "Cadmium";
  vatom_mass[i] *= 112.41;
  vatom_volume[i] = 22.0408;
  vatom_valence_std[i] = 12;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = 4.05;
  vatom_miedema_nws[i] = 1.24;
  vatom_miedema_Vm[i] = 5.5;
  vatom_miedema_gamma_s[i] = 780;
  vatom_miedema_BVm[i] = 6.10;
  vatom_radius[i] = 0.150;
  vatom_radius_covalent[i] = 1.44;
  vatom_electronegativity[i] = 1.69;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 1.36;
  vatom_pearson_coefficient[i] = 0.000271603; // Cd
  // p-electron systems
  i++;
  vatom_symbol[i] = "In";
  vatom_name[i] = "Indium";
  vatom_mass[i] *= 114.82;
  vatom_volume[i] = 27.5233;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.90;
  vatom_miedema_nws[i] = 1.17;
  vatom_miedema_Vm[i] = 6.3;
  vatom_miedema_gamma_s[i] = 690;
  vatom_miedema_BVm[i] = 6.4;
  vatom_radius[i] = 0.157;
  vatom_radius_covalent[i] = 1.42;
  vatom_electronegativity[i] = 1.78;
  vatom_crystal[i] = "fct";
  vatom_pettifor_scale[i] = 1.60;
  vatom_pearson_coefficient[i] = 1.24494E-05; // In
  i++;
  vatom_symbol[i] = "Sn";
  vatom_name[i] = "Tin";
  vatom_mass[i] *= 118.69;
  vatom_volume[i] = 27.5555;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 4.15;
  vatom_miedema_nws[i] = 1.24;
  vatom_miedema_Vm[i] = 6.4;
  vatom_miedema_gamma_s[i] = 710;
  vatom_miedema_BVm[i] = 8.8;
  vatom_radius[i] = 0.158;
  vatom_radius_covalent[i] = 1.39;
  vatom_electronegativity[i] = 1.96;
  vatom_crystal[i] = "bct";
  vatom_pettifor_scale[i] = 1.84;
  vatom_pearson_coefficient[i] = 0.000334085; // Sn
  i++;
  vatom_symbol[i] = "Sb";
  vatom_name[i] = "Antimony";
  vatom_mass[i] *= 121.75;
  vatom_volume[i] = 27.1823;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 4.40;
  vatom_miedema_nws[i] = 1.26;
  vatom_miedema_Vm[i] = 6.6;
  vatom_miedema_gamma_s[i] = 680;
  vatom_miedema_BVm[i] = 7.0;
  vatom_radius[i] = 0.161;
  vatom_radius_covalent[i] = 1.39;
  vatom_electronegativity[i] = 2.05;
  vatom_crystal[i] = "rhl";
  vatom_pettifor_scale[i] = 2.08;
  vatom_pearson_coefficient[i] = 6.60751E-05; // Sb
  i++;
  vatom_symbol[i] = "Te";
  vatom_name[i] = "Tellurium";
  vatom_mass[i] *= 127.6;
  vatom_volume[i] = 28.1993;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 4.72;
  vatom_miedema_nws[i] = 1.31;
  vatom_miedema_Vm[i] = 6.439;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.143;
  vatom_radius_covalent[i] = 1.38;
  vatom_electronegativity[i] = 2.10;
  vatom_crystal[i] = "hex";
  vatom_pettifor_scale[i] = 2.32;
  vatom_pearson_coefficient[i] = 0.000283934; // Te Table 27 of JX
  i++;
  vatom_symbol[i] = "I";
  vatom_name[i] = "Iodine";
  vatom_mass[i] *= 126.9045;
  vatom_volume[i] = 34.9784;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = 5.33;
  vatom_miedema_nws[i] = 0.17;
  vatom_miedema_Vm[i] = 8.72;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.136;
  vatom_radius_covalent[i] = 1.39;
  vatom_electronegativity[i] = 2.66;
  vatom_crystal[i] = "orc";
  vatom_pettifor_scale[i] = 2.56;
  vatom_pearson_coefficient[i] = 0.0; // I interpolation phi_star, nws, Vm,
  i++;
  vatom_symbol[i] = "Xe";
  vatom_name[i] = "Xenon";
  vatom_mass[i] *= 131.3;
  vatom_volume[i] = -1.0000;
  vatom_valence_std[i] = 0;
  vatom_valence_iupac[i] = 8;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.218;
  vatom_radius_covalent[i] = 1.40;
  vatom_electronegativity[i] = 2.60;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.000267781; // Xe JX CHANGED VALENCE

  // ROW6
  // s-electron systems
  i++;
  vatom_symbol[i] = "Cs";
  vatom_name[i] = "Cesium";
  vatom_mass[i] *= 132.9054;
  vatom_volume[i] = 117.281;
  vatom_valence_std[i] = 1;
  vatom_valence_iupac[i] = 1;
  vatom_miedema_phi_star[i] = 1.95;
  vatom_miedema_nws[i] = 0.55;
  vatom_miedema_Vm[i] = 16.8;
  vatom_miedema_gamma_s[i] = 95;
  vatom_miedema_BVm[i] = 1.4;
  vatom_radius[i] = 0.265;
  vatom_radius_covalent[i] = 2.44;
  vatom_electronegativity[i] = 0.79;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.25;
  vatom_pearson_coefficient[i] = 0.0; // Cs
  i++;
  vatom_symbol[i] = "Ba";
  vatom_name[i] = "Barium";
  vatom_mass[i] *= 137.33;
  vatom_volume[i] = 62.6649;
  vatom_valence_std[i] = 2;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = 2.32;
  vatom_miedema_nws[i] = 0.81;
  vatom_miedema_Vm[i] = 11.3;
  vatom_miedema_gamma_s[i] = 370;
  vatom_miedema_BVm[i] = 3.9;
  vatom_radius[i] = 0.217;
  vatom_radius_covalent[i] = 2.15;
  vatom_electronegativity[i] = 0.89;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.50;
  vatom_pearson_coefficient[i] = 6.23705E-05; // Ba
  // d-electron systems: transition metals
  i++;
  vatom_symbol[i] = "La";
  vatom_name[i] = "Lanthanium";
  vatom_mass[i] *= 138.9055;
  vatom_volume[i] = 36.8495;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.05; /*3.17*/
  vatom_miedema_nws[i] = 1.09; /*1.18*/
  vatom_miedema_Vm[i] = 8.0; /*7.98*/
  vatom_miedema_gamma_s[i] = 900;
  vatom_miedema_BVm[i] = 5.5;
  vatom_radius[i] = 0.187;
  vatom_radius_covalent[i] = 2.07;
  vatom_electronegativity[i] = 1.10;
  vatom_crystal[i] = "hex";
  vatom_pettifor_scale[i] = 0.7480;
  vatom_pearson_coefficient[i] = 4.65323E-08; // La
  // lantanidies
  i++;
  vatom_symbol[i] = "Ce";
  vatom_name[i] = "Cerium";
  vatom_mass[i] *= 140.12;
  vatom_volume[i] = 26.4729;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.18;
  vatom_miedema_nws[i] = 1.19;
  vatom_miedema_Vm[i] = 7.76;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.182;
  vatom_radius_covalent[i] = 2.04;
  vatom_electronegativity[i] = 1.12;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0.7460;
  vatom_pearson_coefficient[i] = 2.24956E-05; // Ce pettifor linear interpolation // miedema from Alonso-March.
  i++;
  vatom_symbol[i] = "Pr";
  vatom_name[i] = "Praseodymium";
  vatom_mass[i] *= 140.9077;
  vatom_volume[i] = 36.4987;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.19;
  vatom_miedema_nws[i] = 1.20;
  vatom_miedema_Vm[i] = 7.56;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.183;
  vatom_radius_covalent[i] = 2.03;
  vatom_electronegativity[i] = 1.13;
  vatom_crystal[i] = "hex";
  vatom_pettifor_scale[i] = 0.7440;
  vatom_pearson_coefficient[i] = 0.0; // Pr pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Nd";
  vatom_name[i] = "Neodymium";
  vatom_mass[i] *= 144.24;
  vatom_volume[i] = 29.6719;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.19;
  vatom_miedema_nws[i] = 1.20;
  vatom_miedema_Vm[i] = 7.51;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.182;
  vatom_radius_covalent[i] = 2.01;
  vatom_electronegativity[i] = 1.14;
  vatom_crystal[i] = "hex";
  vatom_pettifor_scale[i] = 0.7420;
  vatom_pearson_coefficient[i] = 0.000231599; // Nd pettifor linear interpolation JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Pm";
  vatom_name[i] = "Promethium";
  vatom_mass[i] *= 146.92;
  vatom_volume[i] = 34.6133;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.19;
  vatom_miedema_nws[i] = 1.21;
  vatom_miedema_Vm[i] = 7.43;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 1.99;
  vatom_electronegativity[i] = 1.13;
  vatom_crystal[i] = "hex";
  vatom_pettifor_scale[i] = 0.7400;
  vatom_pearson_coefficient[i] = 0.0; // Pm pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Sm";
  vatom_name[i] = "Samarium";
  vatom_mass[i] *= 150.4;
  vatom_volume[i] = 33.9484;
  vatom_valence_std[i] = 8;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.20;
  vatom_miedema_nws[i] = 1.21;
  vatom_miedema_Vm[i] = 7.37;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.181;
  vatom_radius_covalent[i] = 1.98;
  vatom_electronegativity[i] = 1.17;
  vatom_crystal[i] = "rhl";
  vatom_pettifor_scale[i] = 0.7380;
  vatom_pearson_coefficient[i] = 0.000334686; // Sm pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Eu";
  vatom_name[i] = "Europium";
  vatom_mass[i] *= 151.96;
  vatom_volume[i] = 43.1719;
  vatom_valence_std[i] = 9;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.20;
  vatom_miedema_nws[i] = 1.21;
  vatom_miedema_Vm[i] = 7.36;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.204;
  vatom_radius_covalent[i] = 1.98;
  vatom_electronegativity[i] = 1.20;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.7360;
  vatom_pearson_coefficient[i] = 4.32857E-05; // Eu pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Gd";
  vatom_name[i] = "Gadolinium";
  vatom_mass[i] *= 157.25;
  vatom_volume[i] = 32.5777;
  vatom_valence_std[i] = 10;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.20;
  vatom_miedema_nws[i] = 1.21;
  vatom_miedema_Vm[i] = 7.34;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.180;
  vatom_radius_covalent[i] = 1.96;
  vatom_electronegativity[i] = 1.20;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.7340;
  vatom_pearson_coefficient[i] = 0.000127674; // Gd pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Tb";
  vatom_name[i] = "Terbium";
  vatom_mass[i] *= 158.9254;
  vatom_volume[i] = 32.0200;
  vatom_valence_std[i] = 11;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.21;
  vatom_miedema_nws[i] = 1.22;
  vatom_miedema_Vm[i] = 7.20;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.177;
  vatom_radius_covalent[i] = 1.94;
  vatom_electronegativity[i] = 1.10;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.7320;
  vatom_pearson_coefficient[i] = 0.0; // Tb pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Dy";
  vatom_name[i] = "Dysprosium";
  vatom_mass[i] *= 162.5;
  vatom_volume[i] = 31.5096;
  vatom_valence_std[i] = 12;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.21;
  vatom_miedema_nws[i] = 1.22;
  vatom_miedema_Vm[i] = 7.12;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.177;
  vatom_radius_covalent[i] = 1.92;
  vatom_electronegativity[i] = 1.22;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.7300;
  vatom_pearson_coefficient[i] = 5.20771E-05; // Dy pettifor linear interpolation JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Ho";
  vatom_name[i] = "Holmium";
  vatom_mass[i] *= 164.9304;
  vatom_volume[i] = 31.0155;
  vatom_valence_std[i] = 13;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.22;
  vatom_miedema_nws[i] = 1.22;
  vatom_miedema_Vm[i] = 7.06;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.176;
  vatom_radius_covalent[i] = 1.92;
  vatom_electronegativity[i] = 1.23;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.7280;
  vatom_pearson_coefficient[i] = 2.96961E-32; // Ho pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Er";
  vatom_name[i] = "Erbium";
  vatom_mass[i] *= 167.26;
  vatom_volume[i] = 30.5431;
  vatom_valence_std[i] = 14;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.22;
  vatom_miedema_nws[i] = 1.23;
  vatom_miedema_Vm[i] = 6.98;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.175;
  vatom_radius_covalent[i] = 1.89;
  vatom_electronegativity[i] = 1.24;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.7260;
  vatom_pearson_coefficient[i] = 7.24618E-05; // Er pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Tm";
  vatom_name[i] = "Thulium";
  vatom_mass[i] *= 168.9342;
  vatom_volume[i] = 30.0016;
  vatom_valence_std[i] = 15;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.22;
  vatom_miedema_nws[i] = 1.23;
  vatom_miedema_Vm[i] = 6.90;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.174;
  vatom_radius_covalent[i] = 1.90;
  vatom_electronegativity[i] = 1.25;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.7240;
  vatom_pearson_coefficient[i] = 0.0; // Tm pettifor linear interpolation JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Yb";
  vatom_name[i] = "Ytterbium";
  vatom_mass[i] *= 173.04;
  vatom_volume[i] = 39.4395;
  vatom_valence_std[i] = 16;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.22;
  vatom_miedema_nws[i] = 1.23;
  vatom_miedema_Vm[i] = 6.86;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.193;
  vatom_radius_covalent[i] = 1.87;
  vatom_electronegativity[i] = 1.10;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0.7220;
  vatom_pearson_coefficient[i] = 8.54557E-05; // Yb pettifor linear interpolation
  i++;
  vatom_symbol[i] = "Lu";
  vatom_name[i] = "Lutetium";
  vatom_mass[i] *= 174.967;
  vatom_volume[i] = 29.3515;
  vatom_valence_std[i] = 17;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.22;
  vatom_miedema_nws[i] = 1.24;
  vatom_miedema_Vm[i] = 6.81;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.173;
  vatom_radius_covalent[i] = 1.87;
  vatom_electronegativity[i] = 1.27;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.7200;
  vatom_pearson_coefficient[i] = 8.27273E-07; // Lu
  // d-electron systems: transition metalsnnn";
  i++;
  vatom_symbol[i] = "Hf";
  vatom_name[i] = "Hafnium";
  vatom_mass[i] *= 178.49;
  vatom_volume[i] = 22.0408;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.55;
  vatom_miedema_nws[i] = 1.43;
  vatom_miedema_Vm[i] = 5.6;
  vatom_miedema_gamma_s[i] = 2200;
  vatom_miedema_BVm[i] = 15.0;
  vatom_radius[i] = 0.159;
  vatom_radius_covalent[i] = 1.75;
  vatom_electronegativity[i] = 1.30;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.775;
  vatom_pearson_coefficient[i] = 5.25384E-05; // Hf
  i++;
  vatom_symbol[i] = "Ta";
  vatom_name[i] = "Tantalum";
  vatom_mass[i] *= 180.9479;
  vatom_volume[i] = 18.1100;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 4.05;
  vatom_miedema_nws[i] = 1.63;
  vatom_miedema_Vm[i] = 4.9;
  vatom_miedema_gamma_s[i] = 3050;
  vatom_miedema_BVm[i] = 22.0;
  vatom_radius[i] = 0.147;
  vatom_radius_covalent[i] = 1.70;
  vatom_electronegativity[i] = 1.50;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.83;
  vatom_pearson_coefficient[i] = 3.66845E-09; // Ta
  i++;
  vatom_symbol[i] = "W";
  vatom_name[i] = "Tungsten";
  vatom_mass[i] *= 183.85;
  vatom_volume[i] = 15.9387;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 4.80;
  vatom_miedema_nws[i] = 1.81;
  vatom_miedema_Vm[i] = 4.5;
  vatom_miedema_gamma_s[i] = 3300;
  vatom_miedema_BVm[i] = 31.0;
  vatom_radius[i] = 0.137;
  vatom_radius_covalent[i] = 1.62;
  vatom_electronegativity[i] = 2.36;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0.885;
  vatom_pearson_coefficient[i] = 6.96679E-05; // W
  i++;
  vatom_symbol[i] = "Re";
  vatom_name[i] = "Rhenium";
  vatom_mass[i] *= 186.2;
  vatom_volume[i] = 14.8941;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = 5.40;
  vatom_miedema_nws[i] = 1.86;
  vatom_miedema_Vm[i] = 4.3;
  vatom_miedema_gamma_s[i] = 3650;
  vatom_miedema_BVm[i] = 33.0;
  vatom_radius[i] = 0.138;
  vatom_radius_covalent[i] = 1.51;
  vatom_electronegativity[i] = 1.90;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.94;
  vatom_pearson_coefficient[i] = 2.70849E-05; // Re
  i++;
  vatom_symbol[i] = "Os";
  vatom_name[i] = "Osmium";
  vatom_mass[i] *= 190.2;
  vatom_volume[i] = 14.2403;
  vatom_valence_std[i] = 8;
  vatom_valence_iupac[i] = 8;
  vatom_miedema_phi_star[i] = 5.40;
  vatom_miedema_nws[i] = 1.85;
  vatom_miedema_Vm[i] = 4.2;
  vatom_miedema_gamma_s[i] = 3500;
  vatom_miedema_BVm[i] = 35.0;
  vatom_radius[i] = 0.135;
  vatom_radius_covalent[i] = 1.44;
  vatom_electronegativity[i] = 2.20;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 0.995;
  vatom_pearson_coefficient[i] = 7.45234E-05; // Os JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Ir";
  vatom_name[i] = "Iridium";
  vatom_mass[i] *= 192.22;
  vatom_volume[i] = 14.5561;
  vatom_valence_std[i] = 9;
  vatom_valence_iupac[i] = 8;
  vatom_miedema_phi_star[i] = 5.55;
  vatom_miedema_nws[i] = 1.83;
  vatom_miedema_Vm[i] = 4.2;
  vatom_miedema_gamma_s[i] = 3100;
  vatom_miedema_BVm[i] = 25.0;
  vatom_radius[i] = 0.135;
  vatom_radius_covalent[i] = 1.41;
  vatom_electronegativity[i] = 2.20;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.05;
  vatom_pearson_coefficient[i] = 2.53787E-05; // Ir JX CHANGED VALENCE
  i++;
  vatom_symbol[i] = "Pt";
  vatom_name[i] = "Platinum";
  vatom_mass[i] *= 195.09;
  vatom_volume[i] = 15.7298;
  vatom_valence_std[i] = 10;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = 5.65;
  vatom_miedema_nws[i] = 1.78;
  vatom_miedema_Vm[i] = 4.4;
  vatom_miedema_gamma_s[i] = 2550;
  vatom_miedema_BVm[i] = 18.0;
  vatom_radius[i] = 0.138;
  vatom_radius_covalent[i] = 1.36;
  vatom_electronegativity[i] = 2.28;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.105;
  vatom_pearson_coefficient[i] = 3.39206E-05; // Pt
  i++;
  vatom_symbol[i] = "Au";
  vatom_name[i] = "Gold";
  vatom_mass[i] *= 196.9665;
  vatom_volume[i] = 18.1904;
  vatom_valence_std[i] = 11;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 5.15;
  vatom_miedema_nws[i] = 1.57;
  vatom_miedema_Vm[i] = 4.7;
  vatom_miedema_gamma_s[i] = 1550;
  vatom_miedema_BVm[i] = 18.0;
  vatom_radius[i] = 0.144;
  vatom_radius_covalent[i] = 1.36;
  vatom_electronegativity[i] = 2.54;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.16;
  vatom_xray_scatt[i] = 74.99;
  vatom_pearson_coefficient[i] = 2.08217E-32; // Au
  i++;
  vatom_symbol[i] = "Hg";
  vatom_name[i] = "Mercury";
  vatom_mass[i] *= 200.59;
  vatom_volume[i] = 29.7156;
  vatom_valence_std[i] = 12;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 4.20;
  vatom_miedema_nws[i] = 1.24;
  vatom_miedema_Vm[i] = 5.8;
  vatom_miedema_gamma_s[i] = 610;
  vatom_miedema_BVm[i] = 4.0;
  vatom_radius[i] = 0.150;
  vatom_radius_covalent[i] = 1.32;
  vatom_electronegativity[i] = 2.00;
  vatom_crystal[i] = "rhl";
  vatom_pettifor_scale[i] = 1.32;
  vatom_pearson_coefficient[i] = 6.52519E-05; // Hg
  // p-electron systems
  i++;
  vatom_symbol[i] = "Tl";
  vatom_name[i] = "Thallium";
  vatom_mass[i] *= 204.37;
  vatom_volume[i] = 31.0721;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = 3.90;
  vatom_miedema_nws[i] = 1.12;
  vatom_miedema_Vm[i] = 6.6;
  vatom_miedema_gamma_s[i] = 610;
  vatom_miedema_BVm[i] = 6.2;
  vatom_radius[i] = 0.171;
  vatom_radius_covalent[i] = 1.45;
  vatom_electronegativity[i] = 1.62;
  vatom_crystal[i] = "hcp";
  vatom_pettifor_scale[i] = 1.56;
  vatom_pearson_coefficient[i] = 1.99659E-05; // Tl electronegativity  2.04=>1.62
  i++;
  vatom_symbol[i] = "Pb";
  vatom_name[i] = "Lead";
  vatom_mass[i] *= 207.2;
  vatom_volume[i] = 31.6649;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 4.10;
  vatom_miedema_nws[i] = 1.15;
  vatom_miedema_Vm[i] = 6.9;
  vatom_miedema_gamma_s[i] = 610;
  vatom_miedema_BVm[i] = 7.9;
  vatom_radius[i] = 0.175;
  vatom_radius_covalent[i] = 1.46;
  vatom_electronegativity[i] = 2.33;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 1.80;
  vatom_pearson_coefficient[i] = 1.94378E-05; // Pb
  i++;
  vatom_symbol[i] = "Bi";
  vatom_name[i] = "Bismuth";
  vatom_mass[i] *= 208.9804;
  vatom_volume[i] = 31.5691;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = 4.15;
  vatom_miedema_nws[i] = 1.16;
  vatom_miedema_Vm[i] = 7.2;
  vatom_miedema_gamma_s[i] = 550;
  vatom_miedema_BVm[i] = 6.7;
  vatom_radius[i] = 0.182;
  vatom_radius_covalent[i] = 1.48;
  vatom_electronegativity[i] = 2.02;
  vatom_crystal[i] = "rhl";
  vatom_pettifor_scale[i] = 2.04;
  vatom_pearson_coefficient[i] = 0.0; // Bi
  i++;
  vatom_symbol[i] = "Po";
  vatom_name[i] = "Polonium";
  vatom_mass[i] *= 209.98;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.140;
  vatom_radius_covalent[i] = 1.40;
  vatom_electronegativity[i] = 2.00;
  vatom_crystal[i] = "sc";
  vatom_pettifor_scale[i] = 2.28;
  vatom_pearson_coefficient[i] = 0.0; // Po
  i++;
  vatom_symbol[i] = "At";
  vatom_name[i] = "Astatine";
  vatom_mass[i] *= 210;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 1.50;
  vatom_electronegativity[i] = 2.20;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 2.52;
  vatom_pearson_coefficient[i] = 0.0; // At
  i++;
  vatom_symbol[i] = "Rn";
  vatom_name[i] = "Radon";
  vatom_mass[i] *= 222;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 0;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 1.50;
  vatom_electronegativity[i] = 2.2;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Rn

  // ROW7
  // s-electron systems
  i++;
  vatom_symbol[i] = "Fr";
  vatom_name[i] = "Francium";
  vatom_mass[i] *= 223.02;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 1;
  vatom_valence_iupac[i] = 1;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 2.60;
  vatom_electronegativity[i] = 0.70;
  vatom_crystal[i] = "bcc";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Fr
  i++;
  vatom_symbol[i] = "Ra";
  vatom_name[i] = "Radium";
  vatom_mass[i] *= 226.0254;
  vatom_volume[i] = -1.0000;
  vatom_valence_std[i] = 2;
  vatom_valence_iupac[i] = 2;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 2.21;
  vatom_electronegativity[i] = 0.89;
  vatom_crystal[i] = "bct";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Ra
  // d-electron systems: transition metals
  i++;
  vatom_symbol[i] = "Ac";
  vatom_name[i] = "Actinium";
  vatom_mass[i] *= 227.03;
  vatom_volume[i] = 45.2437;
  vatom_valence_std[i] = 3;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 2.15;
  vatom_electronegativity[i] = 1.10;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Ac
  // actinidies
  i++;
  vatom_symbol[i] = "Th";
  vatom_name[i] = "Thorium";
  vatom_mass[i] *= 232.0381;
  vatom_volume[i] = 31.9586;
  vatom_valence_std[i] = 4;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = 3.30;
  vatom_miedema_nws[i] = 1.28;
  vatom_miedema_Vm[i] = 7.3;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.180;
  vatom_radius_covalent[i] = 2.06;
  vatom_electronegativity[i] = 1.30;
  vatom_crystal[i] = "fcc";
  vatom_pettifor_scale[i] = 0;
  vatom_xray_scatt[i] = 86.64;
  vatom_pearson_coefficient[i] = 0.0; // Th
  i++;
  vatom_symbol[i] = "Pa";
  vatom_name[i] = "Protoactinium";
  vatom_mass[i] *= 231.04;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 5;
  vatom_valence_iupac[i] = 5;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 2.00;
  vatom_electronegativity[i] = 1.50;
  vatom_crystal[i] = "bct";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Pa
  i++;
  vatom_symbol[i] = "U";
  vatom_name[i] = "Uranium";
  vatom_mass[i] *= 238.03;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 6;
  vatom_valence_iupac[i] = 6;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = 0.138;
  vatom_radius_covalent[i] = 1.96;
  vatom_electronegativity[i] = 1.38;
  vatom_crystal[i] = "orc";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 1.15611E-06; // U
  i++;
  vatom_symbol[i] = "Np";
  vatom_name[i] = "Neptunium";
  vatom_mass[i] *= 237.05;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 7;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 1.90;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Np
  i++;
  vatom_symbol[i] = "Pu";
  vatom_name[i] = "Plutonium";
  vatom_mass[i] *= 244.06;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 8;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 1.87;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Pu
  i++;
  vatom_symbol[i] = "Am";
  vatom_name[i] = "Americium";
  vatom_mass[i] *= 243.06;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 9;
  vatom_valence_iupac[i] = 7;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 1.80;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Am
  i++;
  vatom_symbol[i] = "Cm";
  vatom_name[i] = "Curium";
  vatom_mass[i] *= 247.07;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 10;
  vatom_valence_iupac[i] = 8;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = 1.69;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Cm

  i++;
  vatom_symbol[i] = "Bk";
  vatom_name[i] = "Berkelium";
  vatom_mass[i] *= 247.07;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 11;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = NNN;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Bk
  i++;
  vatom_symbol[i] = "Cf";
  vatom_name[i] = "Californium";
  vatom_mass[i] *= 251.08;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 12;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = NNN;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Cf
  i++;
  vatom_symbol[i] = "Es";
  vatom_name[i] = "Einsteinium";
  vatom_mass[i] *= 252.08;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 13;
  vatom_valence_iupac[i] = 4;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = NNN;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Es
  i++;
  vatom_symbol[i] = "Fm";
  vatom_name[i] = "Fermium";
  vatom_mass[i] *= 257.1;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 14;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = NNN;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Fm
  i++;
  vatom_symbol[i] = "Md";
  vatom_name[i] = "Mendelevium";
  vatom_mass[i] *= 258.1;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 15;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = NNN;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Md
  i++;
  vatom_symbol[i] = "No";
  vatom_name[i] = "Nobelium";
  vatom_mass[i] *= 259.1;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 16;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = NNN;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // No
  i++;
  vatom_symbol[i] = "Lr";
  vatom_name[i] = "Lawrencium";
  vatom_mass[i] *= 262.11;
  vatom_volume[i] = NNN;
  vatom_valence_std[i] = 17;
  vatom_valence_iupac[i] = 3;
  vatom_miedema_phi_star[i] = NNN;
  vatom_miedema_nws[i] = NNN;
  vatom_miedema_Vm[i] = NNN;
  vatom_miedema_gamma_s[i] = NNN;
  vatom_miedema_BVm[i] = NNN;
  vatom_radius[i] = NNN;
  vatom_radius_covalent[i] = NNN;
  vatom_electronegativity[i] = NNN;
  vatom_crystal[i] = "nnn";
  vatom_pettifor_scale[i] = 0;
  vatom_pearson_coefficient[i] = 0.0; // Lr

  //   int valence_WSETYAWAN[]={9999,1,0,1,2,3,4,-3,-2,-1,0,1,2,3,4,-3,-2,-1,0,1,2,3,4,3,3,2,3,3,2,1,2,3,4,-3,-2,-1,0,1,2,3,4,2,4,7,4,3,2,1,2,3,4,-3,-2,-1,0,1,2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,3,4,5,3,4,4,4,2,1,2,3,4,3,0,0,0,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
  //   for(uint j=1;j<93;j++) {
  //   int test=vatom_valence_std[j];
  //   if(j<=9 && test>=5) {test=test-8;} // first row
  //   if(j<=17 && test>=5) {test=test-8;} // second row
  ////   if(test!=valence_WSETYAWAN[j])
  // cerr << vatom_symbol[j] << " " << vatom_valence_std[j] << " " <<  test  << " " << valence_WSETYAWAN[j];
  // if(test!=valence_WSETYAWAN[j]) cerr << "   ****   ";
  // cerr << endl;
  // }
  // throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize():","Throw for debugging purposes.",_GENERIC_ERROR_);

  //   for(int i=1;i<90;i++) {
  //     cerr << i << " " <<  vatom_symbol[i] << " " << vatom_pettifor_scale[i] << endl;
  //   }
  // cerr << i << " " << vatom_symbol[i] << endl; throw aurostd::xerror(__AFLOW_FILE__,XPID+"atoms_initialize","Throw for debugging purposes.",_GENERIC_ERROR_);

  // aconvasp stuff
  // All data collected from the NIST online tables:
  // http://physics.nist.gov/PhysRefData/FFast/html/form.html
  // All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  // These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).
  vatom_xray_scatt[1 + 2] = 3.00145E+00; // Li
  vatom_xray_scatt[1 + 14] = 1.53133E+01; // P
  vatom_xray_scatt[1 + 24] = 2.43589E+01; // Mn
  vatom_xray_scatt[1 + 25] = 2.46830E+01; // Fe

  // All data collected from the online tables:
  // http://www-cxro.lbl.gov/optical_constants/pert_form.html
  // All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
  vatom_xray_scatt[1 + 0] = 1.000; // H
  vatom_xray_scatt[1 + 1] = 2.000; // He
  vatom_xray_scatt[1 + 2] = 3.001; // Li
  vatom_xray_scatt[1 + 6] = 6.019; // C
  vatom_xray_scatt[1 + 7] = 8.052; // O
  vatom_xray_scatt[1 + 13] = 14.43; // P
  vatom_xray_scatt[1 + 14] = 15.30; // P
  vatom_xray_scatt[1 + 20] = 21.34; // Sc
  vatom_xray_scatt[1 + 21] = 22.24; // Ti
  vatom_xray_scatt[1 + 23] = 23.84; // Cr
  vatom_xray_scatt[1 + 24] = 24.46; // Mn
  vatom_xray_scatt[1 + 25] = 24.85; // Fe
  vatom_xray_scatt[1 + 26] = 24.59; // Co
  vatom_xray_scatt[1 + 27] = 25.02; // Ni
  vatom_xray_scatt[1 + 28] = 27.03; // Cu
  vatom_xray_scatt[1 + 29] = 28.44; // Zn
  vatom_xray_scatt[1 + 46] = 47.18; // Ag
  vatom_xray_scatt[1 + 78] = 74.99; // Au
  vatom_xray_scatt[1 + 89] = 86.64; // Th

  // Atomic masses
  // All indices are the atomic number shifted back by one.
  // All masses are in kilograms

  // not useful anymore all the masses are declared
  //  vatom_mass=vector<double> (NUM_ELEMENTS,0.0);
  //  for(i=0;i<NUM_ELEMENTS;i++) {
  //    vatom_mass[i]=(double)(2*i)*AMU2KILOGRAM;
  //  }
  //  vatom_mass[1+0]=1.0079*AMU2KILOGRAM; //H
  //  vatom_mass[1+7]=15.9994*AMU2KILOGRAM; //O
  //  vatom_mass[1+24]=54.93805*AMU2KILOGRAM; //Mn

  // finish and copy
}

// **************************************************************************
// Function GetAtomNumber
// **************************************************************************
uint GetAtomNumber(const string& symbol) {
  for (uint iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return iat;
    }
  }
  cerr << "GetAtomNumber symbol not found symbol=" << symbol << "  vatom_name.size()=" << vatom_name.size() << endl;
  return 0; // no symbol found
}

// **************************************************************************
// Function GetAtomName
// **************************************************************************
std::string GetAtomName(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_name.at(iat);
    }
  }
  return symbol;
}

std::string GetAtomName(const uint& atnum) {
  if (atnum >= vatom_name.size() || atnum <= 0) {
    cerr << "GetAtomName out of boundary  atnum=" << atnum << "  vatom_name.size()=" << vatom_name.size() << endl;
    return "not found";
  }
  return vatom_name.at(atnum);
}

// **************************************************************************
// Function GetAtomSymbol
// **************************************************************************
std::string GetAtomSymbol(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_symbol.at(iat)) {
      return vatom_symbol.at(iat);
    }
  }
  return symbol;
}

std::string GetAtomSymbol(const uint& atnum) {
  if (atnum >= vatom_symbol.size() || atnum <= 0) {
    cerr << "GetAtomSymbol out of boundary  atnum=" << atnum << "  vatom_symbol.size()=" << vatom_symbol.size() << endl;
    return "not found";
  }
  return vatom_symbol.at(atnum);
}

// **************************************************************************
// Function GetAtomMass
// **************************************************************************
double GetAtomMass(const string& _symbol, bool clean) {
  // CO20181128
  string symbol = _symbol; // CO20181128
  if (clean) {
    symbol = KBIN::VASP_PseudoPotential_CleanName(symbol);
  } // CO20181128
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_mass.at(iat);
    }
  }
  return (double) NNN;
}

double GetAtomMass(const uint& atnum) {
  if (atnum >= vatom_mass.size() || atnum <= 0) {
    cerr << "GetAtomMass out of boundary  atnum=" << atnum << "  vatom_mass.size()=" << vatom_mass.size() << endl;
    return (double) NNN;
  }
  return vatom_mass.at(atnum);
}

// **************************************************************************
// Function GetAtomComptonCrossSection
// **************************************************************************
double GetAtomComptonCrossSection(const string& symbol) {
  //  cout << "symbol=" << symbol << endl;
  //  cout << "GetAtomNumber(symbol)=" << GetAtomNumber(symbol) << endl;
  return GetAtomComptonCrossSection(GetAtomNumber(symbol));
}

double GetAtomComptonCrossSection(const uint& atnum) {
  // sigma_c at 662 KeV, compton cross section in barn (1 barn = 1e-28 m^2) [ref: Ortiz, Comp Mat Sci 44, 1042 (2009)]
  const double sigma_c[] = {0.29,  0.57,  0.86,  1.15,  1.43,  1.72,  2.00,  2.29,  2.58,  2.86,  3.15,  3.43,  3.72,  4.00,  4.29,  4.57,  4.86,  5.14,  5.43,  5.71,  5.99,  6.28,  6.56,
                            6.84,  7.13,  7.41,  7.69,  7.98,  8.26,  8.54,  8.82,  9.11,  9.39,  9.67,  9.95,  10.24, 10.52, 10.79, 11.08, 11.36, 11.64, 11.92, 12.20, 12.48, 12.76, 13.04,
                            13.32, 13.60, 13.88, 14.16, 14.44, 14.72, 15.00, 15.27, 15.55, 15.83, 16.11, 16.38, 16.66, 16.94, 17.22, 17.50, 17.77, 18.05, 18.33, 18.60, 18.88, 19.16, 19.44,
                            19.71, 19.98, 20.26, 20.54, 20.81, 21.08, 21.36, 21.64, 21.91, 22.19, 22.46, 22.73, 23.01, 23.28, 23.55, 23.82, 24.10, 24.38, 24.64, 24.92, -1.000};
  vector<double> vsigma_c;
  vsigma_c.push_back(NNN);
  for (uint i = 0; sigma_c[i] > 0; i++) {
    vsigma_c.push_back(sigma_c[i]);
  }
  if (atnum >= vsigma_c.size() || atnum <= 0) {
    cerr << "GetAtomComptonCrossSection out of boundary  atnum=" << atnum << "  vsigma_c.size()=" << vsigma_c.size() << endl;
    return (double) NNN;
  }
  return vsigma_c.at(atnum);
}

// **************************************************************************
// Function GetAtomPhotoelectricCrossSection
// **************************************************************************
double GetAtomPhotoelectricCrossSection(const string& symbol) {
  //  cout << "symbol=" << symbol << endl;
  //  cout << "GetAtomNumber(symbol)=" << GetAtomNumber(symbol) << endl;
  return GetAtomPhotoelectricCrossSection(GetAtomNumber(symbol));
}

double GetAtomPhotoelectricCrossSection(const uint& atnum) {
  // sigma_pe photoelectric cross section in barn (1 barn = 1e-28 m^2) %ref: Ortiz, Comp Mat Sci 44, 1042 (2009)
  const double sigma_pe[] = {8.79e-09, 2.70e-07, 2.62e-06, 1.33e-05, 4.59e-05, 1.22e-04, 2.75e-04, 5.51e-04, 9.85e-04, 1.69e-03, 2.62e-03, 3.97e-03, 6.10e-03, 8.46e-03, 1.18e-02, 1.69e-02, 2.19e-02, 2.87e-02,
                             3.68e-02, 4.91e-02, 5.94e-02, 7.74e-02, 9.47e-02, 1.16e-01, 1.40e-01, 1.59e-01, 1.94e-01, 2.21e-01, 2.65e-01, 3.10e-01, 3.56e-01, 4.14e-01, 4.68e-01, 5.49e-01, 6.24e-01, 7.09e-01,
                             7.98e-01, 9.16e-01, 1.02,     1.10,     1.27,     1.49,     1.57,     1.74,     1.92,     2.15,     2.36,     2.60,     2.87,     3.06,     3.30,     3.66,     4.06,     4.31,
                             4.68,     5.03,     5.47,     5.86,     6.41,     6.95,     7.43,     7.96,     8.56,     9.02,     9.78,     1.05e+01, 1.11e+01, 1.19e+01, 1.27e+01, 1.37e+01, 1.44e+01, 1.48e+01,
                             1.57e+01, 1.65e+01, 1.78e+01, 1.94e+01, 2.06e+01, 2.11e+01, 2.23e+01, 2.36e+01, 2.56e+01, 2.65e+01, 2.80e+01, 3.00e+01, 3.15e+01, 3.33e+01, 3.47e+01, 3.65e+01, 3.82e+01, -1.000};
  vector<double> vsigma_pe;
  vsigma_pe.push_back(NNN);
  for (uint i = 0; sigma_pe[i] > 0; i++) {
    vsigma_pe.push_back(sigma_pe[i]);
  }
  if (atnum >= vsigma_pe.size() || atnum <= 0) {
    cerr << "GetAtomPhotoelectricCrossSection out of boundary  atnum=" << atnum << "  vsigma_pe.size()=" << vsigma_pe.size() << endl;
    return (double) NNN;
  }
  return vsigma_pe.at(atnum);
}

// **************************************************************************
// Function GetAtomVolume
// **************************************************************************
double GetAtomVolume(const string& _symbol, bool clean) {
  // CO20181128
  string symbol = _symbol; // CO20181128
  if (clean) {
    symbol = KBIN::VASP_PseudoPotential_CleanName(symbol);
  } // CO20181128
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_volume.at(iat);
    }
  }
  return (double) NNN;
}

double GetAtomVolume(const uint& atnum) {
  if (atnum >= vatom_volume.size() || atnum <= 0) {
    cerr << "GetAtomVolume out of boundary  atnum=" << atnum << "  vatom_volume.size()=" << vatom_volume.size() << endl;
    return (double) NNN;
  }
  return vatom_volume.at(atnum);
}

// **************************************************************************
// Function GetAtomValenceIupac
// **************************************************************************
int GetAtomValenceIupac(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_valence_iupac.at(iat);
    }
  }
  return (int) NNN;
}

int GetAtomValenceIupac(const uint& atnum) {
  if (atnum >= vatom_valence_iupac.size() || atnum <= 0) {
    cerr << "GetAtomValenceIupac out of boundary  atnum=" << atnum << "  vatom_valence_iupac.size()=" << vatom_valence_iupac.size() << endl;
    return (int) NNN;
  }
  return vatom_valence_iupac.at(atnum);
}

// **************************************************************************
// Function GetAtomValenceStd
// **************************************************************************
int GetAtomValenceStd(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_valence_std.at(iat);
    }
  }
  return (int) NNN;
}

int GetAtomValenceStd(const uint& atnum) {
  if (atnum >= vatom_valence_std.size() || atnum <= 0) {
    cerr << "GetAtomValenceStd out of boundary  atnum=" << atnum << "  vatom_valence_std.size()=" << vatom_valence_std.size() << endl;
    return (int) NNN;
  }
  return vatom_valence_std.at(atnum);
}

// **************************************************************************
// Function GetAtomRadius
// **************************************************************************
double GetAtomRadius(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_radius.at(iat);
    }
  }
  return (double) NNN;
}

double GetAtomRadius(const uint& atnum) {
  if (atnum >= vatom_radius.size() || atnum <= 0) {
    cerr << "GetAtomRadius out of boundary  atnum=" << atnum << "  vatom_radius.size()=" << vatom_radius.size() << endl;
    return (double) NNN;
  }
  return vatom_radius.at(atnum);
}

// DX+CO20170904 START
//  **************************************************************************
//  Function GetAtomRadiusCovalent
//  **************************************************************************
double GetAtomRadiusCovalent(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_radius_covalent.at(iat);
    }
  }
  return (double) NNN;
}

double GetAtomRadiusCovalent(const uint& atnum) {
  if (atnum >= vatom_radius_covalent.size() || atnum <= 0) {
    cerr << "GetAtomRadiusCovalent out of boundary  atnum=" << atnum << "  vatom_radius_covalent.size()=" << vatom_radius_covalent.size() << endl;
    return (double) NNN;
  }
  return vatom_radius_covalent.at(atnum);
}

// DX+CO20170904 END

// **************************************************************************
// Function GetAtomElectronegativity
// **************************************************************************
double GetAtomElectronegativity(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_electronegativity.at(iat);
    }
  }
  return (double) NNN;
}

double GetAtomElectronegativity(const uint& atnum) {
  if (atnum >= vatom_electronegativity.size() || atnum <= 0) {
    cerr << "GetAtomElectronegativity out of boundary  atnum=" << atnum << "  vatom_electronegativity.size()=" << vatom_electronegativity.size() << endl;
    return (double) NNN;
  }
  return vatom_electronegativity.at(atnum);
}

// **************************************************************************
// Function GetAtomCrystal
// **************************************************************************
string GetAtomCrystal(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_crystal.at(iat);
    }
  }
  return (string) "nnn";
}

string GetAtomCrystal(const uint& atnum) {
  if (atnum >= vatom_crystal.size() || atnum <= 0) {
    cerr << "GetAtomCrystal out of boundary  atnum=" << atnum << "  vatom_crystal.size()=" << vatom_crystal.size() << endl;
    return (string) "nnn";
  }
  return vatom_crystal.at(atnum);
}

// **************************************************************************
// Function GetAtomPettiforScale
// **************************************************************************
double GetAtomPettiforScale(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_pettifor_scale.at(iat);
    }
  }
  return 0.0;
}

double GetAtomPettiforScale(const uint& atnum) {
  if (atnum >= vatom_crystal.size() || atnum <= 0) {
    cerr << "GetAtomPettiforScale out of boundary  atnum=" << atnum << "  vatom_crystal.size()=" << vatom_crystal.size() << endl;
    return 0.0;
  }
  return vatom_pettifor_scale.at(atnum);
}

bool GetAtomPettiforScale(const vector<string>& vsymbol, vector<double>& vvalue) {
  vvalue.clear(); // delete
  for (size_t i = 0; i < vsymbol.size(); i++) {
    vvalue.push_back(GetAtomPettiforScale(vsymbol[i]));
  }
  return true;
}

bool GetAtomPettiforScale(const vector<uint>& vatnum, vector<double>& vvalue) {
  vvalue.clear(); // delete
  for (size_t i = 0; i < vatnum.size(); i++) {
    vvalue.push_back(GetAtomPettiforScale(vatnum[i]));
  }
  return true;
}

bool GetAtomPettiforScale(const vector<string>& vsymbol, xvector<double>& vvalue) {
  if (vvalue.rows != (int) vsymbol.size()) {
    return false; // nothing to be ordered
  }
  if (vvalue.lrows != 1) {
    return false; // start from 1
  }
  for (size_t i = 0; i < vsymbol.size(); i++) {
    vvalue[i + 1] = GetAtomPettiforScale(vsymbol[i]);
  }
  return true;
}

bool GetAtomPettiforScale(const vector<uint>& vatnum, xvector<double>& vvalue) {
  if (vvalue.rows != (int) vatnum.size()) {
    return false; // nothing to be ordered
  }
  if (vvalue.lrows != 1) {
    return false; // start from 1
  }
  for (size_t i = 0; i < vatnum.size(); i++) {
    vvalue[i + 1] = GetAtomPettiforScale(vatnum[i]);
  }
  return true;
}

bool SortAtomsPettiforScale(vector<string>& vsymbol, xvector<int>& vorder, xvector<double>& vvalue) {
  if (vorder.rows != (int) vsymbol.size()) {
    return false; // nothing to be ordered
  }
  if (vvalue.rows != (int) vsymbol.size()) {
    return false; // nothing to be ordered
  }
  if (vorder.lrows != 1) {
    return false; // start from 1 .. and contains order from 1
  }
  if (vvalue.lrows != 1) {
    return false; // start from 1
  }
  // build
  for (size_t i = 1; i <= vsymbol.size(); i++) {
    vvalue[i] = GetAtomPettiforScale(vsymbol.at(i - 1));
    vorder[i] = i;
  }
  aurostd::sort2(vsymbol.size(), vvalue, vorder);
  vector<string> vsymbol_tmp(vsymbol);
  vsymbol.clear();
  for (size_t i = 1; i <= vsymbol_tmp.size(); i++) {
    vsymbol.push_back(vsymbol_tmp.at(vorder[i] - 1));
  }
  return true;
}

bool SortAtomsPettiforScale(vector<string>& vsymbol, vector<int>& vorder, vector<double>& vvalue) {
  //   vorder.clear();vvalue.clear();
  //   xvector<int> xvorder(vsymbol.size());
  //   xvector<double> xvvalue(vsymbol.size());
  //   SortAtomsPettiforScale(vsymbol,xvorder,xvvalue);
  //   for(size_t i=0;i<vsymbol.size();i++) {
  //     vvalue.push_back(xvvalue(i+1));
  //     vorder.push_back(xvorder(i+1));
  //   }
  vector<string> vsymbol_tmp(vsymbol);
  vorder.clear();
  vvalue.clear();
  for (size_t i = 0; i < vsymbol.size(); i++) {
    vvalue.push_back(GetAtomPettiforScale(vsymbol[i]));
    vorder.push_back(i);
  }
  aurostd::sort(vvalue, vorder);
  vsymbol.clear();
  for (size_t i = 0; i < vsymbol_tmp.size(); i++) {
    vsymbol.push_back(vsymbol_tmp.at(vorder.at(i)));
  }
  return true;
}

bool SortAtomsPettiforScale(vector<string>& vsymbol, vector<int>& vorder) {
  vector<double> vvalue;
  SortAtomsPettiforScale(vsymbol, vorder, vvalue);
  return true;
}

bool SortAtomsPettiforScale(vector<string>& vsymbol, vector<double>& vvalue) {
  vector<int> vorder;
  SortAtomsPettiforScale(vsymbol, vorder, vvalue);
  return true;
}

bool SortAtomsPettiforScale(vector<string>& vsymbol) {
  vector<double> vvalue;
  vector<int> vorder;
  SortAtomsPettiforScale(vsymbol, vorder, vvalue);
  return true;
}

// **************************************************************************
// Function GetAtomXrayScatt
// **************************************************************************
double GetAtomXrayScatt(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol.at(iat) || symbol == vatom_name.at(iat)) {
      return vatom_xray_scatt.at(iat);
    }
  }
  return (double) NNN;
}

double GetAtomXrayScatt(const uint& atnum) {
  if (atnum >= vatom_xray_scatt.size() || atnum <= 0) {
    cerr << "GetAtomXrayScatt out of boundary  atnum=" << atnum << "  vatom_xray_scatt.size()=" << vatom_xray_scatt.size() << endl;
    return (double) NNN;
  }
  return vatom_xray_scatt.at(atnum);
}

// DX20181220 - get group of atoms - START
//  **************************************************************************
//  Function GetGroupOfAtoms
//  **************************************************************************
vector<string> GetGroupOfAtoms(string& group_name) {
  vector<string> element_list;
  if (aurostd::tolower(group_name) == "metals") {
    element_list.emplace_back("Li");
    element_list.emplace_back("Be");
    element_list.emplace_back("Na");
    element_list.emplace_back("Mg");
    element_list.emplace_back("Al");
    element_list.emplace_back("K");
    element_list.emplace_back("Ca");
    element_list.emplace_back("Sc");
    element_list.emplace_back("Ti");
    element_list.emplace_back("V");
    element_list.emplace_back("Cr");
    element_list.emplace_back("Mn");
    element_list.emplace_back("Fe");
    element_list.emplace_back("Co");
    element_list.emplace_back("Ni");
    element_list.emplace_back("Cu");
    element_list.emplace_back("Zn");
    element_list.emplace_back("Ga");
    element_list.emplace_back("Rb");
    element_list.emplace_back("Sr");
    element_list.emplace_back("Y");
    element_list.emplace_back("Zr");
    element_list.emplace_back("Nb");
    element_list.emplace_back("Mo");
    element_list.emplace_back("Tc");
    element_list.emplace_back("Ru");
    element_list.emplace_back("Rh");
    element_list.emplace_back("Pd");
    element_list.emplace_back("Ag");
    element_list.emplace_back("Cd");
    element_list.emplace_back("In");
    element_list.emplace_back("Sn");
    element_list.emplace_back("Cs");
    element_list.emplace_back("Ba");
    element_list.emplace_back("La");
    element_list.emplace_back("Ce");
    element_list.emplace_back("Pr");
    element_list.emplace_back("Nd");
    element_list.emplace_back("Pm");
    element_list.emplace_back("Sm");
    element_list.emplace_back("Eu");
    element_list.emplace_back("Gd");
    element_list.emplace_back("Tb");
    element_list.emplace_back("Dy");
    element_list.emplace_back("Ho");
    element_list.emplace_back("Er");
    element_list.emplace_back("Tm");
    element_list.emplace_back("Yb");
    element_list.emplace_back("Lu");
    element_list.emplace_back("Hf");
    element_list.emplace_back("Ta");
    element_list.emplace_back("W");
    element_list.emplace_back("Re");
    element_list.emplace_back("Os");
    element_list.emplace_back("Ir");
    element_list.emplace_back("Pt");
    element_list.emplace_back("Au");
    element_list.emplace_back("Hg");
    element_list.emplace_back("Tl");
    element_list.emplace_back("Pb");
    element_list.emplace_back("Bi");
    element_list.emplace_back("Po");
  }
  if (aurostd::tolower(group_name) == "alkali") {
    element_list.emplace_back("Li");
    element_list.emplace_back("Na");
    element_list.emplace_back("K");
    element_list.emplace_back("Rb");
    element_list.emplace_back("Cs");
  }
  if (aurostd::tolower(group_name) == "alkaline") {
    element_list.emplace_back("Be");
    element_list.emplace_back("Mg");
    element_list.emplace_back("Ca");
    element_list.emplace_back("Sr");
    element_list.emplace_back("Ba");
  }
  if (aurostd::tolower(group_name) == "transition_metals" || aurostd::tolower(group_name) == "transition-metals") {
    element_list.emplace_back("Sc");
    element_list.emplace_back("Ti");
    element_list.emplace_back("V");
    element_list.emplace_back("Cr");
    element_list.emplace_back("Mn");
    element_list.emplace_back("Fe");
    element_list.emplace_back("Co");
    element_list.emplace_back("Ni");
    element_list.emplace_back("Cu");
    element_list.emplace_back("Zn");
    element_list.emplace_back("Y");
    element_list.emplace_back("Zr");
    element_list.emplace_back("Nb");
    element_list.emplace_back("Mo");
    element_list.emplace_back("Tc");
    element_list.emplace_back("Ru");
    element_list.emplace_back("Rh");
    element_list.emplace_back("Pd");
    element_list.emplace_back("Ag");
    element_list.emplace_back("Cd");
    element_list.emplace_back("Hf");
    element_list.emplace_back("Ta");
    element_list.emplace_back("W");
    element_list.emplace_back("Re");
    element_list.emplace_back("Os");
    element_list.emplace_back("Ir");
    element_list.emplace_back("Pt");
    element_list.emplace_back("Au");
    element_list.emplace_back("Hg");
    element_list.emplace_back("Tl");
  }
  if (aurostd::tolower(group_name) == "lanthanides") {
    element_list.emplace_back("La");
    element_list.emplace_back("Ce");
    element_list.emplace_back("Pr");
    element_list.emplace_back("Nd");
    element_list.emplace_back("Pm");
    element_list.emplace_back("Sm");
    element_list.emplace_back("Eu");
    element_list.emplace_back("Gd");
    element_list.emplace_back("Tb");
    element_list.emplace_back("Dy");
    element_list.emplace_back("Ho");
    element_list.emplace_back("Er");
    element_list.emplace_back("Tm");
    element_list.emplace_back("Yb");
    element_list.emplace_back("Lu");
  }
  if (aurostd::tolower(group_name) == "nonmetals" || aurostd::tolower(group_name) == "non-metals" || aurostd::tolower(group_name) == "non_metals") {
    cerr << "in non metals" << endl;
    element_list.emplace_back("He");
    element_list.emplace_back("B");
    element_list.emplace_back("C");
    element_list.emplace_back("N");
    element_list.emplace_back("O");
    element_list.emplace_back("F");
    element_list.emplace_back("Ne");
    element_list.emplace_back("Si");
    element_list.emplace_back("P");
    element_list.emplace_back("S");
    element_list.emplace_back("Cl");
    element_list.emplace_back("Ar");
    element_list.emplace_back("As");
    element_list.emplace_back("Se");
    element_list.emplace_back("Br");
    element_list.emplace_back("Kr");
    element_list.emplace_back("Te");
    element_list.emplace_back("I");
    element_list.emplace_back("Xe");
  }
  print(element_list);
  return element_list;
}

// DX20181220 - get group of atoms - END

// **************************************************************************
// Function GetPearsonCoefficient
// **************************************************************************
double GetPearsonCoefficient(const string& symbol) {
  for (int iat = 0; iat < NUM_ELEMENTS; iat++) {
    if (symbol == vatom_symbol[iat] || symbol == vatom_name[iat]) {
      return GetPearsonCoefficient(iat);
    }
  }
  // If not found throw xerror
  const string message = symbol + " is not a valid element name or symbol.";
  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
}

double GetPearsonCoefficient(const int& iat) {
  return vatom_pearson_coefficient.at(iat);
}

// **************************************************************************
// Function GetCompoundAttenuationLength
// **************************************************************************
double GetCompoundAttenuationLength(const vector<string>& species, const vector<double>& composition, const double& density) {
  // density in g/cm^3, return in cm
  if (species.size() != composition.size()) {
    stringstream message; // CO20190629
    message << "species.size()[" << species.size() << "]!=composition.size()[" << composition.size() << "]";
    // CO20190629
    throw aurostd::xerror(__AFLOW_FILE__, "GetCompoundAttenuationLength():", message, _INDEX_MISMATCH_); // CO20190629
  }
  // cout << "Density=" << density << "<br>" << endl;
  double numerator = 0.0;
  double denominator = 0.0;
  for (size_t i = 0; i < species.size(); i++) {
    numerator += composition.at(i) * GetAtomMass(species.at(i)) * 1000; // from Kg to grams
    denominator += density * (composition.at(i) * (GetAtomComptonCrossSection(species[i]) + GetAtomPhotoelectricCrossSection(species[i])) * 1e-24);
  }
  //  cout << "numerator=" << numerator << endl;
  //  cout << "denominator=" << denominator << endl;
  return numerator / denominator; // in cm
}

double GetCompoundAttenuationLength(const deque<string>& _species, const deque<int>& _composition, const double& density) {
  // density in g/cm^3, return in cm
  vector<double> composition;
  for (size_t i = 0; i < _composition.size(); i++) {
    composition.push_back(double(_composition[i]));
  }
  vector<string> species;
  for (size_t i = 0; i < _species.size(); i++) {
    species.emplace_back(_species[i]);
  }
  return GetCompoundAttenuationLength(species, composition, density);
}

// **************************************************************************
// Function XATOM_AlphabetizationSpecies & XATOM_AlphabetizationCompound
// **************************************************************************
string XATOM_AlphabetizationSpecies(const string& speciesA, const string& speciesB) {
  string system;
  if (speciesA <= speciesB) {
    system = speciesA + speciesB;
  } else {
    system = speciesB + speciesA;
  }
  return system;
}

string XATOM_AlphabetizationSpecies(const vector<string>& vspecies_in) {
  vector<string> vspecies(vspecies_in);
  std::sort(vspecies.begin(), vspecies.end());
  return aurostd::joinWDelimiter(vspecies, "");
}

string XATOM_AlphabetizationSpecies(const vector<string>& vspecies_in, const vector<double>& vnumbers_in) {
  stringstream system;
  vector<string> vspecies(vspecies_in);
  vector<double> vnumbers(vnumbers_in);
  aurostd::sort(vspecies, vnumbers);
  for (size_t i = 0; i < vspecies.size(); i++) {
    system << vspecies[i] << vnumbers[i];
  }
  return system.str();
}

void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies, vector<double>& vnumbers) {
  vspecies.clear();
  vnumbers.clear();
  KBIN::VASP_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system), vspecies, vnumbers);
  stringstream systemstream;
  aurostd::sort(vspecies, vnumbers);
  for (size_t i = 0; i < vspecies.size(); i++) {
    systemstream << vspecies[i];
  }
  system = systemstream.str();
}

void XATOM_AlphabetizationCompound(string& system, vector<string>& vspecies, vector<double>& vnumbers) {
  vspecies.clear();
  vnumbers.clear();
  KBIN::VASP_SplitAlloySpecies(KBIN::VASP_PseudoPotential_CleanName(system), vspecies, vnumbers);
  stringstream systemstream;
  aurostd::sort(vspecies, vnumbers);
  for (size_t i = 0; i < vspecies.size(); i++) {
    systemstream << vspecies[i] << vnumbers.at(i);
  }
  system = systemstream.str();
}

void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies) {
  vector<double> vnumbers(vspecies.size());
  XATOM_AlphabetizationSpecies(system, vspecies, vnumbers);
}

void XATOM_AlphabetizationSpecies(string& system) {
  vector<string> vspecies;
  vector<double> vnumbers;
  XATOM_AlphabetizationSpecies(system, vspecies, vnumbers);
}

void XATOM_AlphabetizationCompound(string& system) {
  vector<string> vspecies;
  vector<double> vnumbers;
  XATOM_AlphabetizationCompound(system, vspecies, vnumbers);
}

// **************************************************************************
// Function XATOM_SplitAlloySpecies
// **************************************************************************
uint XATOM_SplitAlloySpecies(const string& alloy_in, vector<string>& speciesX) {
  speciesX = aurostd::getElements(alloy_in); // CO20200624
  return speciesX.size();
}

uint XATOM_SplitAlloySpecies(const string& alloy_in, vector<string>& speciesX, vector<double>& natomsX) {
  speciesX = aurostd::getElements(alloy_in, natomsX);
  return speciesX.size();
}

// **************************************************************************
// Function XATOM_SplitAlloyPseudoPotentials
// **************************************************************************
uint XATOM_SplitAlloyPseudoPotentials(const string& alloy_in, vector<string>& species_ppX) {
  species_ppX = aurostd::getElements(alloy_in, pp_string, false, false, true);
  // CO20200624 - no clean or sort, but do keep_pp
  return species_ppX.size();
}

uint XATOM_SplitAlloyPseudoPotentials(const string& alloy_in, vector<string>& species_ppX, vector<double>& natomsX) {
  species_ppX = aurostd::getElements(alloy_in, natomsX, pp_string, false, false, true);
  // CO20200624 - no clean or sort, but do keep_pp - will return natomX to be all 1's
  return species_ppX.size();
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// _SYM_OP
// look into aflow.h for the definitions

// constructors
_sym_op::_sym_op() {
  free();
}
_sym_op::_sym_op(const _sym_op& b) {
  copy(b);
}

// destructor
_sym_op::~_sym_op() {
  free();
}

void _sym_op::free() {
  Uc.clear();
  Uf.clear(); // clear stuff
  generator.clear(); // clear stuff
  generator_coefficients.clear(); // clear stuff       //DX20171206 - generator coefficients
  SU2_matrix.clear(); // clear stuff	      //DX20180115 - 2x2 complex SU(2) matrix
  su2_coefficients.clear(); // clear stuff       //DX20180115 - su(2) coefficients on Pauli matrices
  angle = 0.0; // clear stuff
  axis.clear(); // clear stuff
  quaternion_vector.clear(); // clear stuff	//GG
  quaternion_matrix.clear(); // clear stuff	//GG
  str_type = ""; // clear stuff
  str_Hermann_Mauguin = ""; // clear stuff
  str_Schoenflies = ""; // clear stuff
  flag_inversion = false; // clear stuff
  is_pgroup = false; // clear stuff
  is_pgroup_xtal = false; // clear stuff
  is_pgroupk_Patterson = false; // clear stuff      //DX20200129
  is_pgroupk = false; // clear stuff
  is_pgroupk_xtal = false; // clear stuff       //DX20171205
  ctau.clear();
  ftau.clear(); // clear stuff
  basis_atoms_map.clear(); // clear stuff
  basis_types_map.clear(); // clear stuff
  basis_map_calculated = false; // clear stuff
  is_fgroup = false; // clear stuff
  ctrasl.clear();
  ftrasl.clear(); // clear stuff
  is_sgroup = false; // clear stuff
  site = 0; // clear stuff       //DX20170803
  is_agroup = false; // clear stuff
}

void _sym_op::copy(const _sym_op& b) {
  Uc = b.Uc;
  Uf = b.Uf;
  generator = b.generator;
  generator_coefficients = b.generator_coefficients; // DX20171206 - generator coefficients
  SU2_matrix = b.SU2_matrix; // DX20180115 - 2x2 complex SU(2) matrix
  su2_coefficients = b.su2_coefficients; // DX20180115 - su(2) coefficients on Pauli matrices
  angle = b.angle;
  axis = b.axis;
  quaternion_vector = b.quaternion_vector; // GG
  quaternion_matrix = b.quaternion_matrix; // GG
  str_type = b.str_type;
  str_Hermann_Mauguin = b.str_Hermann_Mauguin;
  str_Schoenflies = b.str_Schoenflies;
  flag_inversion = b.flag_inversion;
  is_pgroup = b.is_pgroup;
  is_pgroup_xtal = b.is_pgroup_xtal;
  is_pgroupk_Patterson = b.is_pgroupk_Patterson; // DX20200129
  is_pgroupk = b.is_pgroupk;
  is_pgroupk_xtal = b.is_pgroupk_xtal; // DX20171205
  ctau = b.ctau;
  ftau = b.ftau;
  basis_atoms_map.clear();
  for (size_t i = 0; i < b.basis_atoms_map.size(); i++) {
    basis_atoms_map.push_back(b.basis_atoms_map[i]);
  }
  basis_types_map.clear();
  for (size_t i = 0; i < b.basis_types_map.size(); i++) {
    basis_types_map.push_back(b.basis_types_map[i]);
  }
  basis_map_calculated = b.basis_map_calculated;
  is_fgroup = b.is_fgroup;
  ctrasl = b.ctrasl;
  ftrasl = b.ftrasl;
  is_sgroup = b.is_sgroup;
  site = b.site; // DX20170803
  is_agroup = b.is_agroup;
}

const _sym_op& _sym_op::operator=(const _sym_op& b) {
  // operator=
  if (this != &b) {
    free();
    copy(b);
  }
  return *this;
}

ostream& operator<<(ostream& oss, const _sym_op& symop) {
  xmatrix<double> Uexp(3, 3);
  // oss.setf(std::ios::fixed,std::ios::floatfield);
  // oss.precision(10);
  if (symop.is_pgroup == true) {
    oss << " pgroup" << endl;
  }
  if (symop.is_pgroup_xtal == true) {
    oss << " pgroup_xtal" << endl;
  }
  if (symop.is_pgroupk_Patterson == true) {
    oss << " pgroupk_Patterson" << endl; // DX20200129
  }
  if (symop.is_fgroup == true) {
    oss << " fgroup" << endl;
  }
  if (symop.is_sgroup == true) {
    oss << " sgroup" << endl;
  }
  if (symop.is_agroup == true) {
    oss << " agroup" << endl;
  }
  if (symop.is_pgroupk == true) {
    oss << " pgroupk" << endl;
  }
  if (symop.is_pgroupk_xtal == true) {
    oss << " pgroupk_xtal" << endl; // DX20171205
  }
  if (symop.is_agroup == true) {
    oss << " site:" << symop.site << endl; // DX20170803
  }
  oss << " type: " << symop.str_type << endl;
  oss << " Hermann_Mauguin: " << symop.str_Hermann_Mauguin << endl;
  oss << " Schoenflies: " << symop.str_Schoenflies << endl;
  // DX+CO START
  oss << "" << roundoff(symop.Uc, 1e-8) << " Uc " << endl; // CO roundoff for printing
  // oss << "" << symop.Uc << " Uc "<< endl;
  // DX+CO END
  Uexp = exp(symop.generator);
  // if(! a.symop_inversion[k] ) oss << "" << max(aurostd::abs(Uexp-Uc))
  // << " error " << endl; else oss << "" <<  max(aurostd::abs(Uexp+Uc)) << " error " << endl;
  // oss << "" << xint(symop.Uf) << " Uf "<< endl; //CO, WILL ZERO OUT 0.999999999999999999
  // CO START
  const xmatrix<int> Uf_int(symop.Uf.lrows, symop.Uf.lcols, symop.Uf.urows, symop.Uf.ucols);
  for (int i = symop.Uf.lrows; i <= symop.Uf.urows; i++) {
    for (int j = symop.Uf.lcols; j <= symop.Uf.ucols; j++) {
      Uf_int[i][j] = aurostd::nint(symop.Uf[i][j]);
    }
  }
  oss << "" << Uf_int << " Uf " << endl;
  // CO END
  oss << "" << symop.generator << " A=generator U=+-exp(A) [not Uc and -1 if inversion]" << endl;
  // DX20171206
  oss << "" << symop.generator_coefficients << "  so(3) expansion coefficients on Lx, Ly, and Lz basis" << endl;
  // DX20171206
  // DX20180115 - adding SU(2) and su(2); specific xcomplex printing - START
  char buf11_re[80];
  char buf11_im[80];
  char buf12_re[80];
  char buf12_im[80];
  char buf21_re[80];
  char buf21_im[80];
  char buf22_re[80];
  char buf22_im[80];
  const string iobuf = "%11.4le";
  // HE20221102 switching from deprecated sprintf to snprintf (eliminates the chance of buffer overflows)
  snprintf(buf11_re, 80, iobuf.c_str(), symop.SU2_matrix(1, 1).re);
  snprintf(buf11_im, 80, iobuf.c_str(), symop.SU2_matrix(1, 1).im);
  snprintf(buf12_re, 80, iobuf.c_str(), symop.SU2_matrix(1, 2).re);
  snprintf(buf12_im, 80, iobuf.c_str(), symop.SU2_matrix(1, 2).im);
  snprintf(buf21_re, 80, iobuf.c_str(), symop.SU2_matrix(2, 1).re);
  snprintf(buf21_im, 80, iobuf.c_str(), symop.SU2_matrix(2, 1).im);
  snprintf(buf22_re, 80, iobuf.c_str(), symop.SU2_matrix(2, 2).re);
  snprintf(buf22_im, 80, iobuf.c_str(), symop.SU2_matrix(2, 2).im);
  oss << " (" << buf11_re << "," << buf11_im << ") ";
  oss << " (" << buf12_re << "," << buf12_im << ")" << endl;
  oss << " (" << buf21_re << "," << buf21_im << ") ";
  oss << " (" << buf22_re << "," << buf22_im << ")" << "  SU(2) complex matrix [(real,imaginary)]" << endl;
  // DX - formatting issues with xcomplex: oss << " " << symop.SU2_matrix << "  SU(2) complex matrix [(real,imaginary)]" << endl;
  char buf1_re[80];
  char buf1_im[80];
  char buf2_re[80];
  char buf2_im[80];
  char buf3_re[80];
  char buf3_im[80];
  snprintf(buf1_re, 80, iobuf.c_str(), symop.su2_coefficients(1).re);
  snprintf(buf1_im, 80, iobuf.c_str(), symop.su2_coefficients(1).im);
  snprintf(buf2_re, 80, iobuf.c_str(), symop.su2_coefficients(2).re);
  snprintf(buf2_im, 80, iobuf.c_str(), symop.su2_coefficients(2).im);
  snprintf(buf3_re, 80, iobuf.c_str(), symop.su2_coefficients(3).re);
  snprintf(buf3_im, 80, iobuf.c_str(), symop.su2_coefficients(3).im);
  oss << " (" << buf1_re << "," << buf1_im << ") ";
  oss << " (" << buf2_re << "," << buf2_im << ") ";
  oss << " (" << buf3_re << "," << buf3_im << ")" << "  su(2) expansion coefficients on Pauli matrices [(real,imaginary)]" << endl;
  // DX - formatting issues with complex: oss << " " << symop.su2_coefficients << "  su(2) expansion coefficients on Pauli matrices [(real,imaginary)]" << endl;
  // DX20180115 - adding SU(2) and su(2); specific xcomplex printing - END
  oss << " " << symop.angle << "  angle " << endl;
  oss << "" << symop.axis << "  axis " << endl;
  // GG START HERE
  //  Quaternion Output
  //  Roundoff used to print zeros instead of values e-17
  oss << "" << roundoff(symop.quaternion_vector, 1e-8) << " quaternion_vector " << endl;
  oss << "" << roundoff(symop.quaternion_matrix, 1e-8) << " quaternion_matrix " << endl;
  // GG STOP HERE
  oss << " " << symop.flag_inversion << "   inversion " << endl;
  //   oss << " "<< symop.str_type << endl;
  //   oss << " "<< symop.str_Hermann_Mauguin << endl;
  //   oss << " "<< symop.str_Schoenflies << endl;
  if (symop.is_fgroup || symop.is_sgroup) {
    oss << "" << symop.ctau << "    ctau " << endl;
  }
  if (symop.is_fgroup || symop.is_sgroup) {
    oss << "" << symop.ftau << "    ftau " << endl;
  }
  if (symop.is_sgroup) {
    oss << "" << symop.ctrasl << "    ctrasl " << endl;
  }
  if (symop.is_sgroup) {
    oss << "" << symop.ftrasl << "    ftrasl " << endl;
  }

  // DX+CO START
  // if(symop.is_fgroup==true||symop.is_agroup==true) {
  if (symop.basis_map_calculated) {
    // DX+CO END
    oss << " - ";
    for (size_t n = 0; n < symop.basis_atoms_map.size(); n++) {
      oss << symop.basis_atoms_map[n] << " ";
    }
    oss << "   basis_atoms_map" << " ";
    oss << " - ";
    for (size_t n = 0; n < symop.basis_types_map.size(); n++) {
      oss << symop.basis_types_map[n] << " ";
    }
    oss << "   basis_types_map" << " ";
    oss << endl;
  }

  return oss;
}

void _sym_op::setUc(const xmatrix<double>& _Uc, const xmatrix<double>& lattice) {
  // CO20190520
  Uc = _Uc;
  const xmatrix<double> f2c = trasp(lattice);
  const xmatrix<double> c2f = inverse(trasp(lattice));
  Uf = c2f * Uc * f2c;
}

void _sym_op::setUf(const xmatrix<double>& _Uf, const xmatrix<double>& lattice) {
  // CO20190520
  Uf = _Uf;
  const xmatrix<double> f2c = trasp(lattice);
  const xmatrix<double> c2f = inverse(trasp(lattice));
  Uc = f2c * Uf * c2f;
}

// DX201801107 - add _kpoint class - START
//  ***************************************************************************
//  ***************************************************************************
//  _kpoint

// constructors
_kpoint::_kpoint() {
  iomode = IOAFLOW_AUTO;
  klattice.clear(); // clear stuff
  fpos.clear();
  cpos.clear(); // clear stuff
  label = ""; // clear stuff
  is_transformed = false; // clear stuff
}

// destructor
_kpoint::~_kpoint() {
  free();
}

void _kpoint::free() {}

// assignment operator
const _kpoint& _kpoint::operator=(const _kpoint& b) {
  // operator=
  if (this != &b) {
    free();
    iomode = b.iomode;
    fpos = b.fpos;
    cpos = b.cpos;
    label = b.label;
    is_transformed = b.is_transformed;
  }
  return *this;
}

// print kpoint string
string _kpoint::str() const {
  // CO20220611 - negative sign formatting
  //[CO20220611 - adding spacing for minus]string tmp = "   " + aurostd::joinWDelimiter(xvecDouble2vecString(fpos,4,true,1e-4,FIXED_STREAM),"   ") + "   ! " + label;
  stringstream oss;
  oss.precision(4);
  oss << "   ";
  for (int i = fpos.lrows; i <= fpos.urows; i++) {
    if (fpos[i] < 10) {
      oss << " "; // will never happen
    }
    if (!std::signbit(fpos[i])) {
      oss << " ";
    }
    oss << aurostd::utype2string(fpos[i], 4, true, 1e-4, FIXED_STREAM);
    if (i < fpos.urows) {
      oss << "   ";
    }
  }
  oss << "   ! " << label;
  if (is_transformed) {
    oss << "\'";
  } // add prime
  return oss.str();
}

// operator<<
ostream& operator<<(ostream& oss, const _kpoint& kpt) {
  oss << kpt.str();
  return oss;
}

// transform kpoint
void _kpoint::TransformKpoint(const xmatrix<double>& P) {
  fpos = fpos * P;
  klattice = aurostd::inverse(P) * klattice; // i.e., klattice'=Q*klattice
  is_transformed = true;
}

// DX201801107 - add _kpoint class - END

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// wyckoffsite_ITC
// look into aflow.h for the definitions
// constructors

wyckoffsite_ITC::wyckoffsite_ITC() {
  coord.clear();
  index = 0; // DX20200427
  type = "";
  wyckoffSymbol = "";
  letter = ""; // DX20180128 - add Wyckoff letter
  site_symmetry = ""; // DX20180128 - add Wyckoff site symmetry
  multiplicity = 0; // DX20180128 - add Wyckoff multiplicity
  site_occupation = 1.0; // DX20191010 - add site occupation (default: 1.0)
  equations.clear(); // DX20180128 - add Wyckoff multiplicity
  parameter_index = 0; // DX20200513
}

// destructor
wyckoffsite_ITC::~wyckoffsite_ITC() {
  free();
}

// empty free
void wyckoffsite_ITC::free() {}

// operator=
const wyckoffsite_ITC& wyckoffsite_ITC::operator=(const wyckoffsite_ITC& b) {
  // operator=
  if (this != &b) {
    free();
    coord = b.coord;
    index = b.index; // DX20200501
    type = b.type;
    wyckoffSymbol = b.wyckoffSymbol;
    letter = b.letter; // DX20180128 - add Wyckoff letter
    site_symmetry = b.site_symmetry; // DX20180128 - add Wyckoff site symmetry
    multiplicity = b.multiplicity; // DX20180128 - add Wyckoff multiplicity
    site_occupation = b.site_occupation; // DX20191010 - add site occupation
    equations = b.equations; // DX20180128 - add Wyckoff multiplicity
    parameter_index = b.parameter_index; // DX20200513
  }
  return *this;
}

// DX20190130 - add comparison operator so we can sort by Wyckoff letter, then species - START
//  operator<
bool wyckoffsite_ITC::operator<(const wyckoffsite_ITC& b) const {
  // operator<
  if (letter < b.letter) {
    return true;
  } else if (letter == b.letter) {
    if (type < b.type) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

// DX20190130 - add comparison operator so we can sort by Wyckoff letter, then species - END

// copy
wyckoffsite_ITC::wyckoffsite_ITC(const wyckoffsite_ITC& b) {
  free();
  coord = b.coord;
  index = b.index; // DX20200501
  type = b.type;
  wyckoffSymbol = b.wyckoffSymbol;
  letter = b.letter; // DX20180128 - add Wyckoff letter
  site_symmetry = b.site_symmetry; // DX20180128 - add Wyckoff site symmetry
  multiplicity = b.multiplicity; // DX20180128 - add Wyckoff multiplicity
  site_occupation = b.site_occupation; // DX20191010 - add site occupation
  equations = b.equations; // DX20180128 - add Wyckoff multiplicity
  parameter_index = b.parameter_index; // DX20200513
}

// operator <<
ostream& operator<<(ostream& oss, const wyckoffsite_ITC& site) {
  oss << " coord: " << site.coord << endl;
  oss << " index: " << site.index << endl; // DX20200501
  oss << " type: " << site.type << endl;
  oss << " letter: " << site.letter << endl;
  oss << " site_symmetry: " << site.site_symmetry << endl;
  oss << " multiplicity: " << site.multiplicity << endl;
  oss << " wyckoffSymbol: " << site.wyckoffSymbol << endl;
  oss << " site_occupation: " << site.site_occupation << endl; // DX20191010 - add site occupation
  oss << " parameter_index: " << site.parameter_index << endl; // DX20200513
  oss << " equations: " << endl; // DX20191010 - add site occupation
  for (size_t i = 0; i < site.equations.size(); i++) {
    oss << "  " << aurostd::joinWDelimiter(site.equations[i], ",") << endl;
  }
  return oss;
}

// ***************************************************************************
// Implementing JSON SERIALIZABLE HERE
// ***************************************************************************

// Serializable for_atom
aurostd::JSON::object _atom::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_atom_MEMBERS)});
}

_atom _atom::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_atom_MEMBERS)
  return *this;
}

// Serializable for _sym_op
aurostd::JSON::object _sym_op::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_symop_MEMBERS)});
}

_sym_op _sym_op::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_symop_MEMBERS)
  return *this;
}

// Serializable for wyckoffsite_ITC
aurostd::JSON::object wyckoffsite_ITC::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_wyckoffsite_MEMBERS)});
}

wyckoffsite_ITC wyckoffsite_ITC::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_wyckoffsite_MEMBERS)
  return *this;
}

#endif // _AFLOW_XATOM_CPP
