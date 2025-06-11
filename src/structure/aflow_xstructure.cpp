
#include "structure/aflow_xstructure.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <ios>
#include <iostream>
#include <istream>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_automatic_template.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xtensor.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "aflowlib/aflowlib_entry_loader.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_support_types.h"
#include "flow/aflow_xclasses.h"
#include "interfaces/aflow_voro.h"
#include "modules/COMPARE/aflow_compare_structure.h"
#include "modules/HULL/aflow_chull.h" //HE20210408
#include "modules/PROTOTYPES/aflow_anrl.h"
#include "modules/SYM/aflow_symmetry.h"
#include "modules/SYM/aflow_symmetry_spacegroup.h" //DX20180723
#include "structure/aflow_lattice.h"
#include "structure/aflow_xatom.h"

#define _calculate_symmetry_default_sgroup_radius_ 2.0
#define PLATON_MIN_VOLUME_PER_ATOM 6.0   // for symmetry calculation

#define PLATON_TOLERANCE_ANGLE 1.0
#define PLATON_TOLERANCE_D1 0.25
#define PLATON_TOLERANCE_D2 0.25
#define PLATON_TOLERANCE_D3 0.25

#define _EPS_ 0.02

#define DEBUG_LATTICE_DIMENSIONS false // DX20210401

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ios_base;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

// ***************************************************************************
// getAtomIndicesByType() //DX20210322
// ***************************************************************************
vector<uint> getAtomIndicesByType(const xstructure& xstr, int type) {
  // Get the atom indices of a given type from an xstructure

  stringstream message;

  const size_t natoms = xstr.atoms.size();

  vector<uint> indices_atoms_subset;
  for (uint i = 0; i < natoms; i++) {
    if (xstr.atoms[i].type == type) {
      indices_atoms_subset.push_back(i);
    }
  }

  if (indices_atoms_subset.empty()) {
    message << "No atoms found with type = " << type << ". Check structure.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  return indices_atoms_subset;
}

// ***************************************************************************
// getAtomIndicesByName() //DX20210322
// ***************************************************************************
vector<uint> getAtomIndicesByName(const xstructure& xstr, const string& name) {
  // Get the atom indices of a given name/species from an xstructure

  stringstream message;

  const size_t natoms = xstr.atoms.size();
  const string name_clean = KBIN::VASP_PseudoPotential_CleanName(name);

  vector<uint> indices_atoms_subset;
  for (uint i = 0; i < natoms; i++) {
    if (xstr.atoms[i].cleanname == name_clean) {
      indices_atoms_subset.push_back(i);
    }
  }

  if (indices_atoms_subset.empty()) {
    message << "No atoms found with name = " << name << " (note, using name_clean=" << name_clean << "). Check structure.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  return indices_atoms_subset;
}

// ***************************************************************************
// getLeastFrequentAtomTypes() //DX20210322
// ***************************************************************************
vector<uint> getLeastFrequentAtomTypes(const xstructure& xstr) {
  // The least frequent atom type is the atom type with the smallest
  // concentration in the crystal. The atoms of this type are the minimal set
  // that exhibit the crystal periodicity (useful for finding
  // alternative lattices and translations).
  // Returns ALL LFA types (could be more than one)

  vector<uint> lfa_types; // lfa = least frequent atom

  // find minimum type count
  const int type_count_min = aurostd::min(xstr.num_each_type);

  // find all species with this atom count
  for (size_t i = 0; i < xstr.num_each_type.size(); i++) {
    if (xstr.num_each_type[i] == type_count_min) {
      lfa_types.push_back(i);
    }
  }

  if (lfa_types.empty()) {
    throw aurostd::xerror(__AFLOW_FILE__, "getLeastFrequentAtomTypes():", "Least frequent atom type not found. Bad xstructure.", _INPUT_ERROR_);
  }

  return lfa_types;
}

// ***************************************************************************
// getLeastFrequentAtomSpecies() //DX20201230 - moved from XtalFinder
// ***************************************************************************
vector<string> getLeastFrequentAtomSpecies(const xstructure& xstr, bool clean) {
  // The least frequent atom species is the species with the smallest
  // concentration in the crystal. These atoms are the minimal set
  // that exhibit the crystal periodicity (useful for finding
  // alternative lattices and translations).
  // clean: cleans atom name (removes pseudopotential)
  // Returns ALL LFA species (could be more than one)

  vector<string> lfa_species; // lfa = least frequent atom

  // find minimum type count
  const int type_count_min = aurostd::min(xstr.num_each_type);

  // find the first species with this atom count
  for (size_t i = 0; i < xstr.num_each_type.size(); i++) {
    if (xstr.num_each_type[i] == type_count_min) {
      if (clean) {
        lfa_species.push_back(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i]));
      } else {
        lfa_species.push_back(xstr.species[i]);
      }
    }
  }

  if (lfa_species.empty()) {
    throw aurostd::xerror(__AFLOW_FILE__, "getLeastFrequentAtomTypes():", "Least frequent atom species not found. Bad xstructure.", _INPUT_ERROR_);
  }

  return lfa_species;
}

// **************************************************************************
// Function xstructure::GetElements() //DX20200728
// **************************************************************************
vector<string> xstructure::GetElements(bool clean_name, bool fake_names) const { // Made function const //SD20220221

  // Returns the elements in the xstructure

  const bool LDEBUG = (false || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // 1) try xstructure.species
  if (!species.empty()) { // DX20210315
    if (clean_name) {
      vector<string> vspecies;
      for (size_t i = 0; i < species.size(); i++) {
        vspecies.push_back(KBIN::VASP_PseudoPotential_CleanName(species[i]));
      }
      return vspecies;
    } else {
      return aurostd::deque2vector((*this).species);
    }
  }
  // ---------------------------------------------------------------------------
  // 2) try element names (check if first atom name is given)
  else if (!atoms[0].name.empty()) {
    return GetElementsFromAtomNames(clean_name);
  }
  // ---------------------------------------------------------------------------
  // 3) if all are empty, return fake elements (optional) // Does not change xstructure //SD20220221
  else if (atoms[0].name.empty() && fake_names) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " WARNING: Atoms are not labeled" << endl;
    }
    return pflow::getFakeElements(num_each_type.size());
  }

  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "There are no element names in the structure.", _RUNTIME_ERROR_);
}

// **************************************************************************
// Function xstructure::GetElements() //DX20200728
// **************************************************************************
vector<string> xstructure::GetElementsFromAtomNames(bool clean_name) const { // Made function const //SD20220221

  // Extracts the species from the atom names

  vector<string> species_found;
  if (atoms.empty()) {
    return species_found;
  }
  if (!atoms[0].name_is_given) {
    return species_found;
  }

  uint iat = 0;
  string species_tmp;
  for (size_t i = 0; i < num_each_type.size(); i++) {
    species_tmp = atoms[iat].name; // always the first in the species set
    if (species_tmp.empty()) {
      stringstream message;
      message << "Found 1 non-specified species within structure:" << endl;
      message << (*this) << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    for (int j = 0; j < num_each_type[i]; j++) {
      // check all atoms of the same type have the same name
      if (atoms[iat].name != species_tmp) {
        stringstream message;
        message << "The number of each type and atom names do not agree." << endl;
        message << (*this) << endl;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
      iat++;
    }
    if (clean_name) {
      species_found.push_back(KBIN::VASP_PseudoPotential_CleanName(species_tmp));
    } else {
      species_found.push_back(species_tmp);
    }
  }
  return species_found;
}

// **************************************************************************
// Function xstructure::GetReducedComposition() //DX20200728
// **************************************************************************
vector<uint> xstructure::GetReducedComposition(bool numerical_sort) {
  vector<uint> composition;
  for (size_t i = 0; i < num_each_type.size(); i++) {
    composition.push_back((size_t) num_each_type[i]);
  }
  vector<uint> reduced_composition;

  // sort in numerical order (useful for prototypes)
  if (numerical_sort) {
    std::stable_sort(composition.begin(), composition.end());
  }

  // reduce by GCD
  aurostd::reduceByGCD(composition, reduced_composition);

  return reduced_composition;
}

// ***************************************************************************
// Function getCentroidOfStructure() //DX20200728
// ***************************************************************************
xvector<double> getCentroidOfStructure(const xstructure& xstr, bool use_cpos, bool use_atom_mass) {
  return getCentroidOfStructure(xstr.atoms, use_cpos, use_atom_mass);
}

// useful for molecules/non-periodic structures too
xvector<double> getCentroidOfStructure(const deque<_atom>& atoms, bool use_cpos, bool use_atom_mass) {
  // Calculate centroid in a structure
  // This overload is useful for non-periodic structures as well
  // (e.g., needed for symmetry analysis of molecules)
  // Uses aurostd::getCentroid()

  // ---------------------------------------------------------------------------
  // coordinate type
  vector<xvector<double>> coordinates;
  // cpos
  if (use_cpos) {
    for (size_t i = 0; i < atoms.size(); i++) {
      coordinates.push_back(atoms[i].cpos);
    }
  }
  // fpos
  else {
    for (size_t i = 0; i < atoms.size(); i++) {
      coordinates.push_back(atoms[i].fpos);
    }
  }

  // ---------------------------------------------------------------------------
  // coordinate weights
  vector<double> weights;
  // uses mass of atomic species
  if (use_atom_mass) {
    for (size_t i = 0; i < atoms.size(); i++) {
      weights.push_back(atoms[i].mass);
    }
  }
  // consider as geometric points
  else {
    for (size_t i = 0; i < atoms.size(); i++) {
      weights.push_back(1.0);
    }
  }

  return getCentroid(coordinates, weights);
}

// ***************************************************************************
// Function getCentroidOfStructurePBC() //DX20200728
// ***************************************************************************
xvector<double> getCentroidOfStructurePBC(const xstructure& xstr, bool use_cpos, bool use_atom_mass) {
  return getCentroidOfStructurePBC(xstr.atoms, xstr.lattice, use_cpos, use_atom_mass);
}

xvector<double> getCentroidOfStructurePBC(const deque<_atom>& atoms, xmatrix<double> lattice, bool use_cpos, bool use_atom_mass) {
  // Calculate centroid in a structure with periodic boundary conditions
  // Uses aurostd::getCentroidPBC()

  xmatrix<double> cell;

  // ---------------------------------------------------------------------------
  // coordinate type
  vector<xvector<double>> coordinates;
  // cpos
  if (use_cpos) {
    for (size_t i = 0; i < atoms.size(); i++) {
      coordinates.push_back(atoms[i].cpos);
    }
    cell = lattice; // use lattice
  }
  // fpos
  else {
    for (size_t i = 0; i < atoms.size(); i++) {
      coordinates.push_back(atoms[i].fpos);
    }
    cell = aurostd::eye<double>(3, 3); // use unit cube
  }

  // ---------------------------------------------------------------------------
  // coordinate weights
  vector<double> weights;
  // uses mass of atomic species
  if (use_atom_mass) {
    for (size_t i = 0; i < atoms.size(); i++) {
      weights.push_back(atoms[i].mass);
    }
  }
  // consider as geometric points
  else {
    for (size_t i = 0; i < atoms.size(); i++) {
      weights.push_back(1.0);
    }
  }

  return getCentroidPBC(coordinates, weights, cell);
}

// ***************************************************************************
// Reset dims for RadiusSphereLattice() - DX20191122
// ***************************************************************************
void resetLatticeDimensions(
    const xmatrix<double>& lattice, double radius, xvector<int>& dims, vector<xvector<double>>& l1, vector<xvector<double>>& l2, vector<xvector<double>>& l3, vector<int>& a_index, vector<int>& b_index, vector<int>& c_index) {
  // resets the lattice dimensions (dims) based on radius
  // generates lattice vectors (l1,l2,l3) right away = speed increase
  // stores dimension indices (a_index,b_index,c_index)
  // new dims explore order : zeroth cell to max dims = speed increase
  // (can break early if match is found)
  // Useful for finding the closest neihbors or minimum interatomic distances:
  // once we find a neighbor, update/reduce how far we need to search to find
  // a closer neighbor
  // DX create function date: 20190705

  // ---------------------------------------------------------------------------
  // get new dimensions based on radius
  if (radius <= _ZERO_TOL_) {
    dims[1] = 1;
    dims[2] = 1;
    dims[3] = 1;
  } else {
    dims = LatticeDimensionSphere(lattice, radius);
  }
  // ---------------------------------------------------------------------------
  // clear old
  l1.clear();
  l2.clear();
  l3.clear();
  a_index.clear();
  b_index.clear();
  c_index.clear();

  // ---------------------------------------------------------------------------
  // [NEW] - go from zeroth cell out
  // more likely to find match close to origin, why start so far away

  // push back zeroth cell : dims[1]=dims[2]=dims[3]=0
  l1.push_back(0 * lattice(1));
  a_index.push_back(0);
  l2.push_back(0 * lattice(2));
  b_index.push_back(0);
  l3.push_back(0 * lattice(3));
  c_index.push_back(0);

  // push back 1,-1,2,-2,...dims,-dims
  for (int a = 1; a <= dims[1]; a++) {
    l1.push_back(a * lattice(1));
    a_index.push_back(a);
    l1.push_back(-a * lattice(1));
    a_index.push_back(-a);
  }
  for (int b = 1; b <= dims[2]; b++) {
    l2.push_back(b * lattice(2));
    b_index.push_back(b);
    l2.push_back(-b * lattice(2));
    b_index.push_back(-b);
  }
  for (int c = 1; c <= dims[3]; c++) {
    l3.push_back(c * lattice(3));
    c_index.push_back(c);
    l3.push_back(-c * lattice(3));
    c_index.push_back(-c);
  }
}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() - DX20191122
// ***************************************************************************
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice, double& min_dist, uint& frequency, vector<xvector<double>>& coordinates) {
  // determine the minimum coordination shell of the lattice
  // i.e., find the set of closest neighbors to the origin
  // (overload: uses lattice radius)

  // ---------------------------------------------------------------------------
  // determine necessary search radius
  const double radius = RadiusSphereLattice(lattice);

  minimumCoordinationShellLatticeOnly(lattice, min_dist, frequency, coordinates, radius);
}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() - DX20191122
// ***************************************************************************
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice, double& min_dist, uint& frequency, vector<xvector<double>>& coordinates, double radius) {
  // determine the minimum coordination shell of the lattice
  // i.e., find the set of closest neighbors to the origin
  // (overload: instantiates lattice dimension information)

  // ---------------------------------------------------------------------------
  // instantiate lattice vectors
  vector<xvector<double>> l1;
  vector<xvector<double>> l2;
  vector<xvector<double>> l3;
  vector<int> a_index;
  vector<int> b_index;
  vector<int> c_index;
  xvector<int> dims(3);
  dims[1] = dims[2] = dims[3] = 0; // declare/reset
  resetLatticeDimensions(lattice, radius, dims, l1, l2, l3, a_index, b_index, c_index);

  minimumCoordinationShellLatticeOnly(lattice, dims, l1, l2, l3, a_index, b_index, c_index, min_dist, frequency, coordinates, radius);
}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() - DX20191122
// ***************************************************************************
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
                                         xvector<int>& dims,
                                         vector<xvector<double>>& l1,
                                         vector<xvector<double>>& l2,
                                         vector<xvector<double>>& l3,
                                         vector<int>& a_index,
                                         vector<int>& b_index,
                                         vector<int>& c_index,
                                         double& min_dist,
                                         uint& frequency,
                                         vector<xvector<double>>& coordinates,
                                         double radius) {
  // determine the minimum coordination shell environment of the lattice
  // i.e., find the set of closest neighbors to the origin
  // stores l1, l2, l3, a_index, b_index, and c_index for external use
  // optional "radius" as enables more control over search space
  // (and potential speed up, may not need to search as far as the lattice radius)

  xvector<double> tmp;

  // ---------------------------------------------------------------------------
  // reset lattice dimensions
  resetLatticeDimensions(lattice, radius, dims, l1, l2, l3, a_index, b_index, c_index);

  const double relative_tolerance = 10.0; // coordination shell thickness is ten percent of minimum distance

  // ---------------------------------------------------------------------------
  // loop through lattice vectors (stored before-hand in l1,l2,l3)
  for (size_t m = 0; m < l1.size(); m++) {
    const xvector<double> a_component = l1[m]; // DX : coord1-coord2+a*lattice(1)
    for (size_t n = 0; n < l2.size(); n++) {
      const xvector<double> ab_component = a_component + l2[n]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
      for (size_t p = 0; p < l3.size(); p++) {
        if (!(m == 0 && n == 0 && p == 0)) {
          tmp = ab_component + l3[p]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
          const double tmp_mod = aurostd::modulus(tmp);
          // ---------------------------------------------------------------------------
          // if found a new minimum distance and update coordination/frequency and coordinate
          if (tmp_mod < min_dist) {
            // ---------------------------------------------------------------------------
            // if new distance is close to the original it is the same coordination shell (add to coordination)
            // otherwise, reset coordination shell
            // DX - FIXED TOL (bad for undecorated prototypes) - if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
            if (aurostd::isequal(tmp_mod, min_dist, (min_dist / relative_tolerance))) {
              frequency += 1;
              coordinates.push_back(tmp);
            } // tenth of min_dist
            else {
              frequency = 1;
              coordinates.clear();
              coordinates.push_back(tmp);
            } // initialize
            min_dist = tmp_mod;
            // ---------------------------------------------------------------------------
            // diminishing dims: if minimum distance changed, then we may not need to search as far
            // reset loop and search again based on new minimum distance
            if (!(dims[1] == 1 && dims[2] == 1 && dims[3] == 1)) {
              resetLatticeDimensions(lattice, min_dist, dims, l1, l2, l3, a_index, b_index, c_index);
              m = n = p = 0;
              frequency = 0; // reset
              coordinates.clear(); // DX20210222
            }
          }
          // DX - FIXED TOL (bad for undecorated prototypes) - else if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
          else if (aurostd::isequal(tmp_mod, min_dist, (min_dist / relative_tolerance))) {
            frequency += 1;
            coordinates.push_back(tmp);
          } // tenth of min dist
        }
      }
    }
  }
}

// ***************************************************************************
// minimumCoordinationShell() - DX20191122
// ***************************************************************************
void minimumCoordinationShell(const xstructure& xstr, uint center_index, double& min_dist, uint& frequency, vector<xvector<double>>& coordinates) {
  const string type;

  minimumCoordinationShell(xstr, center_index, min_dist, frequency, coordinates, type);
}

// ***************************************************************************
// minimumCoordinationShell() - DX20191122
// ***************************************************************************
void minimumCoordinationShell(const xstructure& xstr, uint center_index, double& min_dist, uint& frequency, vector<xvector<double>>& coordinates, const string& type) {
  // determine the minimum coordination shell environment
  // "type" enables the search of environments by certain elements/types only
  // (e.g., find the neighborhood of oxygen atoms surrounding a magnesium center)

  xvector<double> tmp_coord;

  // ---------------------------------------------------------------------------
  // instantiate lattice vectors
  vector<xvector<double>> l1;
  vector<xvector<double>> l2;
  vector<xvector<double>> l3;
  vector<int> a_index;
  vector<int> b_index;
  vector<int> c_index;
  xvector<int> dims(3);
  dims[1] = dims[2] = dims[3] = 0; // declare/reset
  // resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);

  const double relative_tolerance = 10.0; // coordination shell thickness is ten percent of minimum distance

  for (size_t ii = 0; ii < xstr.atoms.size(); ii++) {
    // ---------------------------------------------------------------------------
    // if atom ii is not environment center, find minimum distance between center atom ii's images
    if (ii != center_index && (xstr.atoms[ii].name == type || type.empty())) { // DX20191105 - added type==""
      const xvector<double> incell_dist = xstr.atoms[center_index].cpos - xstr.atoms[ii].cpos;
      const double incell_mod = aurostd::modulus(incell_dist);
      if (!(dims[1] == 1 && dims[2] == 1 && dims[3] == 1) && incell_mod != 1e9) {
        resetLatticeDimensions(xstr.lattice, incell_mod, dims, l1, l2, l3, a_index, b_index, c_index);
      }
      // DX20180423 - running vector in each loop saves computations; fewer duplicate operations
      for (size_t m = 0; m < l1.size(); m++) {
        const xvector<double> a_component = incell_dist + l1[m]; // DX : coord1-coord2+a*lattice(1)
        for (size_t n = 0; n < l2.size(); n++) {
          const xvector<double> ab_component = a_component + l2[n]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
          for (size_t p = 0; p < l3.size(); p++) {
            tmp_coord = ab_component + l3[p]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
            const double tmp_mod = aurostd::modulus(tmp_coord);
            if (tmp_mod < min_dist) {
              // DX - FIXED TOL (bad for undecorated prototypes) - if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
              if (aurostd::isequal(tmp_mod, min_dist, (min_dist / relative_tolerance))) {
                frequency += 1;
                coordinates.push_back(tmp_coord);
              } // tenth of min_dist
              else {
                frequency = 1;
                coordinates.clear();
                coordinates.push_back(tmp_coord);
              } // initialize
              min_dist = tmp_mod;
            }
            // DX - FIXED TOL (bad for undecorated prototypes) else if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
            else if (aurostd::isequal(tmp_mod, min_dist, (min_dist / relative_tolerance))) {
              frequency += 1;
              coordinates.push_back(tmp_coord);
            } // tenth of min_dist
          }
        }
      }
    }
    // ---------------------------------------------------------------------------
    // if atom is environment center check its images, but only need to search as
    // far as min_dist or lattice_radius (whichever is smaller)
    else if (ii == center_index && (xstr.atoms[ii].name == type || type.empty())) { // DX20191105 - added type==""
      const double lattice_radius = RadiusSphereLattice(xstr.lattice);
      const double search_radius = min(lattice_radius, min_dist);

      // ---------------------------------------------------------------------------
      // use variant that stores the lattice dimension information so it can be
      // updated for the "minimumCoordinationShell" function
      minimumCoordinationShellLatticeOnly(xstr.lattice, dims, l1, l2, l3, a_index, b_index, c_index, min_dist, frequency, coordinates, search_radius);
    }
  }
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// STRUCTURE
// look into aflow.h for the definitions

void xstructure::free() { // DX20191220 - moved all initializations from constuctor into free()
  iomode = IOAFLOW_AUTO; // what else could we do right now...
  title = ""; // DX20191210
  directory = "";
  prototype = "";
  info = "";
  if (title.empty()) {
    title = "NO_TITLE_GIVEN";
  }
  // num_atoms=0;
  // num_types=0;
  scale = 1.0;
  neg_scale = false;
  scale_second = DEFAULT_POCC_SITE_TOL; // DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; //CO20180409
  neg_scale_second = false; // CO20180409
  scale_third.isentry = false; // CO20170803 - site tol
  scale_third.content_double = DEFAULT_POCC_STOICH_TOL; // DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;  //CO20170803 - site tol
  coord_type[0] = coord_type[1] = 0;
  coord_flag = _COORDS_FRACTIONAL_; // _COORDS_FRACTIONAL_ (0) fractional, _COORDS_CARTESIAN_ (1) cartesian.
  isd = false; // !=0 => Selective dynamics, =0 => no selective dynamics.
  lattice.clear();
  lattice(1, 1) = lattice(2, 2) = lattice(3, 3) = 1.0;
  a = b = c = 1.0;
  alpha = beta = gamma = 90.0;
  klattice = ReciprocalLattice(lattice, scale);
  f2c = trasp(lattice);
  c2f = inverse(trasp(lattice));
  symbolic_math_representation_only = false; // DX20180618
  constrained_symmetry_calculation = false; // DX20180618
  symbolic_math_lattice.clear(); // DX20180618
  num_parameters = 0; // number of parameters ANRL 20180618
  num_lattice_parameters = 0; // number of lattice parameters ANRL 20180618
  prototype_parameter_list.clear(); // prototype parameter list ANRL 20180618
  prototype_parameter_values.clear(); // prototype parameters values ANRL 20180618
  origin.clear();
  // TOLERANCES ------------------------
  equiv_fpos_epsilon = _EQUIV_FPOS_EPS_; // standard but you can change
  // NUM_EACH_TYPE ---------------------
  // ClearSpecies(); //CO20180420
  num_each_type.clear(); // CO20180420 - ClearSpecies()
  comp_each_type.clear(); // CO20180420 - ClearSpecies()
  stoich_each_type.clear(); // CO20171025 //CO20180420 - ClearSpecies()
  // SPECIES ---------------------------
  species.clear();
  species_pp.clear();
  species_pp_type.clear();
  species_pp_version.clear();
  species_pp_ZVAL.clear();
  species_pp_vLDAU.clear();
  species_volume.clear();
  species_mass.clear(); // CO20180420 - ClearSpecies()
  is_vasp4_poscar_format = true;
  is_vasp5_poscar_format = false;
  // ATOMS -----------------------------
  atoms.clear();
  // FLAGS -----------------------------
  primitive_calculated = false; // DX20201005
  Niggli_calculated = false;
  Niggli_avoid = false;
  Minkowski_calculated = false;
  Minkowski_avoid = false;
  LatticeReduction_calculated = false;
  LatticeReduction_avoid = false;
  // LATTICE stuff ---------------------
  Standard_Lattice_calculated = false;
  Standard_Lattice_avoid = false;
  Standard_Lattice_primitive = false;
  Standard_Lattice_conventional = false;
  Standard_Lattice_has_failed = false;
  bravais_lattice_type = "";
  bravais_lattice_variation_type = ""; // WSETYAWAN mod
  bravais_lattice_system = "";
  bravais_lattice_lattice_type = "";
  bravais_lattice_lattice_variation_type = ""; // WSETYAWAN mod
  bravais_lattice_lattice_system = "";
  pearson_symbol = "";
  reciprocal_lattice_type = "";
  reciprocal_lattice_variation_type = ""; // WSETYAWAN mod
  volume_changed_original2new = false; // DX20181024
  transform_coordinates_original2new.clear(); // DX20181024
  transform_coordinates_new2original.clear(); // DX20181024
  rotate_lattice_original2new.clear(); // DX20181024
  rotate_lattice_new2original.clear(); // DX20181024
  // reciprocal_conventional_lattice_type="";
  bravais_superlattice_lattice.clear(); // DX20210209
  bravais_superlattice_type = "";
  bravais_superlattice_variation_type = "";
  bravais_superlattice_system = "";
  pearson_symbol_superlattice = "";
  // GENERAL PURPOSE LABEL -------------
  label_uint = 0;
  label_int = 0;
  label_double = 0;
  // ORDER PARAMETER -------------------
  order_parameter_structure = false;
  order_parameter_atoms.clear();
  order_parameter_orbit = 1; // always orbit of itself
  order_parameter_sum = 0;
  // PARTIAL OCCUPATION -------------------
  partial_occupation_flag = false;
  partial_occupation_site_tol = DEFAULT_POCC_SITE_TOL;
  // DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;     // DEFAULT //CO20180409
  partial_occupation_stoich_tol = DEFAULT_POCC_STOICH_TOL;
  // DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;   // DEFAULT //CO20180409
  partial_occupation_HNF = 0;
  partial_occupation_sublattice.clear();
  // FORCES/POSITIONS ------------------
  qm_calculated = false;
  qm_scale = 1.0;
  qm_lattice.clear();
  qm_lattice(1, 1) = qm_lattice(2, 2) = qm_lattice(3, 3) = 1.0;
  qm_klattice = ReciprocalLattice(qm_lattice, qm_scale);
  qm_f2c = trasp(qm_lattice);
  qm_c2f = inverse(trasp(qm_lattice));
  qm_origin.clear();
  qm_atoms.clear();
  qm_forces.clear();
  qm_forces_write = false;
  qm_positions.clear();
  qm_positions_write = false;
  qm_E_cell = 0.0;
  qm_dE_cell = 0.0;
  qm_H_cell = 0.0;
  qm_PV_cell = 0.0;
  qm_mag_cell = 0.0;
  qm_P = 0.0;
  qm_E_atom = 0.0;
  qm_dE_atom = 0.0;
  qm_H_atom = 0.0;
  qm_PV_atom = 0.0;
  qm_mag_atom = 0.0;
  // KPOINTS ---------------------------
  kpoints_k1 = 0;
  kpoints_k2 = 0;
  kpoints_k3 = 0;
  kpoints_s1 = 0;
  kpoints_s2 = 0;
  kpoints_s3 = 0;
  kpoints_kmax = 0;
  kpoints_kppra = 0;
  kpoints_mode = 0;
  kpoints_kscheme = "";
  // DX+CO START
  dist_nn_min = AUROSTD_NAN; // CO
  // SYMMETRY TOLERANCE ----------------------------
  sym_eps = AUROSTD_NAN; // DX
  sym_eps_calculated = false; // DX, this means that it was calculated and set by the symmetry routines
  sym_eps_change_count = 0; // DX20180222 - added tolerance count specific to structure
  sym_eps_no_scan = false; // DX20210331 - added no scan specific to structure
  // DX+CO END
  //  PGROUP ----------------------------
  pgroup.clear(); // just initialize
  pgroup_calculated = false;
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear(); // just initialize
  pgroup_xtal_calculated = false;
  crystal_family = "";
  crystal_system = "";
  point_group_crystal_class = "";
  point_group_Shoenflies = "";
  point_group_Hermann_Mauguin = "";
  point_group_orbifold = "";
  point_group_type = "";
  point_group_order = "";
  point_group_structure = "";
  // PGROUPK_PATTERSON ---------------------------- //DX20200129
  pgroupk_Patterson.clear(); // just initialize
  pgroupk_Patterson_calculated = false;
  // PGROUPK ----------------------------
  pgroupk.clear(); // just initialize
  pgroupk_calculated = false;
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear(); // just initialize //DX20171205 - Added pgroupk_xtal
  pgroupk_xtal_calculated = false; // DX20171205 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear(); // just initialize
  fgroup_calculated = false;
  // SGROUP ----------------------------
  sgroup_radius = -_calculate_symmetry_default_sgroup_radius_; // symmetry not calculated
  sgroup_radius_dims.clear();
  sgroup.clear(); // just initialize
  sgroup_calculated = false;
  // SITE POINT GROUP ------------------
  agroup_calculated = false;
  for (size_t i = 0; i < agroup.size(); i++) {
    agroup[i].clear();
  }
  agroup.clear();
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated = false;
  for (size_t i = 0; i < iatoms.size(); i++) {
    iatoms[i].clear();
  }
  iatoms.clear();
  // SPACE GROUP WITH PLATON/FINDSYM -----------
  spacegroup = "";
  spacegrouplabel = "";
  spacegroupnumber = 0;
  spacegroupnumberoption = 0;
  spacegroupoption = "";
  is_spacegroup_platon = false;
  is_spacegroup_findsym = false;
  is_spacegroup_aflow = false;
  // SPACE GROUP ITC (RHT)
  crystal_system_ITC = ""; // RHT
  point_group_ITC = ""; // RHT
  // DX+CO START
  bravais_label_ITC = 'X'; // RHT
  lattice_label_ITC = 'X'; // RHT
  space_group_ITC = 0; // RHT
  // DX+CO END
  wyckoff_library_entry_ITC = ""; // RHT
  wyccar_ITC.clear(); // RHT
  standard_lattice_ITC.clear(); // RHT
  standard_basis_ITC.clear(); // RHT
  wyckoff_sites_ITC.clear(); // RHT
  wyckoff_symbols_ITC.clear(); // RHT
  setting_ITC = 0; // DX20170830 - SGDATA
  origin_ITC.clear(); // DX20170830 - SGDATA
  general_position_ITC.clear(); // DX20170830 - SGDATA
  // GRID ATOMS ------------------------
  grid_atoms_calculated = false;
  grid_atoms_dimsL.clear();
  grid_atoms_dimsH.clear();
  grid_atoms.clear(); // just initialize
  grid_atoms_number = 0;
  grid_atoms_sc2pcMap.clear(); // CO20171025
  grid_atoms_pc2scMap.clear(); // CO20171025
  // LIJK OBEJCTS ----------------------
  lijk_calculated = false;
  lijk_table.clear();
  lijk_cpos.clear();
  lijk_fpos.clear();
  lijk_dims.clear();
  // NEIGHBOR ------------------------
  // OUTPUT/ERROR ----------------------
  Niggli_has_failed = false;
  Minkowski_has_failed = false;
  LatticeReduction_has_failed = false;
  write_lattice_flag = false;
  write_klattice_flag = false;
  write_inequivalent_flag = false;
  write_DEBUG_flag = false;
  error_flag = false;
  error_string = "";
  // -----------------------------------
}

// Constructors
xstructure::xstructure(const string& structure_title) { // CO20211122
  (*this).initialize(structure_title); // CO20211122
}

// ifstream/istream
xstructure::xstructure(istream& _input, int _iomode) {
  (*this).initialize(_input, _iomode);
}

xstructure::xstructure(ifstream& _input, int _iomode) {
  (*this).initialize(_input, _iomode);
}

xstructure::xstructure(const stringstream& __input, int _iomode) { // DX20210129 - added const
  (*this).initialize(__input, _iomode); // DX20210129
}

xstructure::xstructure(const string& _input, int _iomode) { // CO20211122
  (*this).initialize(_input, _iomode); // CO20211122
}

xstructure::xstructure(const string& url, const string& file, int _iomode) { // CO20211122
  (*this).initialize(url, file, _iomode); // CO20211122
}

void xstructure::copy(const xstructure& bstr) {
  // All the other stuff not set in the constructor
  iomode = bstr.iomode;
  title = bstr.title;
  directory = bstr.directory;
  prototype = bstr.prototype;
  info = bstr.info;
  //  num_types=bstr.num_types;
  //  num_atoms=bstr.num_atoms;
  scale = bstr.scale;
  neg_scale = bstr.neg_scale;
  scale_second = bstr.scale_second; // CO20180409
  neg_scale_second = bstr.neg_scale_second; // CO20180409
  scale_third = bstr.scale_third; // CO20170803 - site tol //CO20180409
  strcpy(coord_type, bstr.coord_type);
  coord_flag = bstr.coord_flag;
  isd = bstr.isd;
  lattice = bstr.lattice;
  a = bstr.a;
  b = bstr.b;
  c = bstr.c;
  alpha = bstr.alpha;
  beta = bstr.beta;
  gamma = bstr.gamma;
  klattice = bstr.klattice;
  origin = bstr.origin;
  f2c = bstr.f2c;
  c2f = bstr.c2f;
  symbolic_math_representation_only = bstr.symbolic_math_representation_only; // DX20180618
  constrained_symmetry_calculation = bstr.constrained_symmetry_calculation; // DX20180618
  symbolic_math_lattice = bstr.symbolic_math_lattice; // DX20180618
  num_parameters = bstr.num_parameters; // number of parameters ANRL 20180618
  num_lattice_parameters = bstr.num_lattice_parameters; // number of lattice parameters ANRL 20180618
  prototype_parameter_list = bstr.prototype_parameter_list; // prototype parameter list ANRL 20180618
  prototype_parameter_values = bstr.prototype_parameter_values; // prototype parameters values ANRL 20180618
  // TOLERANCES ------------------------
  equiv_fpos_epsilon = bstr.equiv_fpos_epsilon;
  // NUM_EACH_TYPE ---------------------
  num_each_type.clear();
  for (size_t i = 0; i < bstr.num_each_type.size(); i++) {
    num_each_type.push_back(bstr.num_each_type[i]);
  }
  comp_each_type.clear();
  for (size_t i = 0; i < bstr.comp_each_type.size(); i++) {
    comp_each_type.push_back(bstr.comp_each_type.at(i));
  }
  stoich_each_type.clear(); // CO20171025
  for (size_t i = 0; i < bstr.stoich_each_type.size(); i++) { // CO20171025
    stoich_each_type.push_back(bstr.stoich_each_type[i]); // CO20171025
  }
  // SPECIES ---------------------------
  species.clear();
  for (size_t i = 0; i < bstr.species.size(); i++) {
    species.push_back(bstr.species[i]);
  }
  species_pp.clear();
  for (size_t i = 0; i < bstr.species_pp.size(); i++) {
    species_pp.push_back(bstr.species_pp[i]);
  }
  species_pp_type.clear();
  for (size_t i = 0; i < bstr.species_pp_type.size(); i++) {
    species_pp_type.push_back(bstr.species_pp_type[i]);
  }
  species_pp_version.clear();
  for (size_t i = 0; i < bstr.species_pp_version.size(); i++) {
    species_pp_version.push_back(bstr.species_pp_version[i]);
  }
  species_pp_ZVAL.clear();
  for (size_t i = 0; i < bstr.species_pp_ZVAL.size(); i++) {
    species_pp_ZVAL.push_back(bstr.species_pp_ZVAL[i]);
  }
  species_pp_vLDAU.clear();
  for (size_t i = 0; i < bstr.species_pp_vLDAU.size(); i++) {
    species_pp_vLDAU.push_back(bstr.species_pp_vLDAU[i]);
  }
  species_volume.clear();
  for (size_t i = 0; i < bstr.species_volume.size(); i++) {
    species_volume.push_back(bstr.species_volume[i]);
  }
  species_mass.clear();
  for (size_t i = 0; i < bstr.species_mass.size(); i++) {
    species_mass.push_back(bstr.species_mass[i]);
  }
  is_vasp4_poscar_format = bstr.is_vasp4_poscar_format;
  is_vasp5_poscar_format = bstr.is_vasp5_poscar_format;
  // FLAGS -----------------------------
  primitive_calculated = bstr.primitive_calculated; // DX20201007
  Niggli_calculated = bstr.Niggli_calculated;
  Niggli_avoid = bstr.Niggli_avoid;
  Minkowski_calculated = bstr.Minkowski_calculated;
  Minkowski_avoid = bstr.Minkowski_avoid;
  LatticeReduction_calculated = bstr.LatticeReduction_calculated;
  LatticeReduction_avoid = bstr.LatticeReduction_avoid;
  // LATTICE stuff ---------------------
  Standard_Lattice_calculated = bstr.Standard_Lattice_calculated;
  Standard_Lattice_avoid = bstr.Standard_Lattice_avoid;
  Standard_Lattice_primitive = bstr.Standard_Lattice_primitive;
  Standard_Lattice_conventional = bstr.Standard_Lattice_conventional;
  Standard_Lattice_has_failed = bstr.Standard_Lattice_has_failed;
  bravais_lattice_type = bstr.bravais_lattice_type;
  bravais_lattice_variation_type = bstr.bravais_lattice_variation_type;
  bravais_lattice_system = bstr.bravais_lattice_system;
  bravais_lattice_lattice_type = bstr.bravais_lattice_lattice_type;
  bravais_lattice_lattice_variation_type = bstr.bravais_lattice_lattice_variation_type;
  bravais_lattice_lattice_system = bstr.bravais_lattice_lattice_system;
  pearson_symbol = bstr.pearson_symbol;
  reciprocal_lattice_type = bstr.reciprocal_lattice_type;
  reciprocal_lattice_variation_type = bstr.reciprocal_lattice_variation_type;
  bravais_superlattice_lattice = bstr.bravais_superlattice_lattice; // DX20210209
  bravais_superlattice_type = bstr.bravais_superlattice_type;
  bravais_superlattice_variation_type = bstr.bravais_superlattice_variation_type;
  bravais_superlattice_system = bstr.bravais_superlattice_system;
  pearson_symbol_superlattice = bstr.pearson_symbol_superlattice;
  volume_changed_original2new = bstr.volume_changed_original2new; // DX20181024
  transform_coordinates_original2new = bstr.transform_coordinates_original2new; // DX20181024
  transform_coordinates_new2original = bstr.transform_coordinates_new2original; // DX20181024
  rotate_lattice_original2new = bstr.rotate_lattice_original2new; // DX20181024
  rotate_lattice_new2original = bstr.rotate_lattice_new2original; // DX20181024
  // ATOMS -----------------------------
  atoms.clear();
  for (size_t i = 0; i < bstr.atoms.size(); i++) {
    atoms.push_back(bstr.atoms[i]);
  }
  // GENERAL PURPOSE LABEL -------------
  label_uint = bstr.label_uint;
  label_int = bstr.label_int;
  label_double = bstr.label_double;
  // ORDER PARAMETER -------------------
  order_parameter_structure = bstr.order_parameter_structure;
  order_parameter_atoms.clear();
  for (size_t i = 0; i < bstr.order_parameter_atoms.size(); i++) {
    order_parameter_atoms.push_back(bstr.order_parameter_atoms[i]);
  }
  order_parameter_orbit = bstr.order_parameter_orbit;
  order_parameter_sum = bstr.order_parameter_sum;
  // PARTIAL OCCUPATION -------------------
  partial_occupation_flag = bstr.partial_occupation_flag;
  partial_occupation_site_tol = bstr.partial_occupation_site_tol; // CO20180409
  partial_occupation_stoich_tol = bstr.partial_occupation_stoich_tol; // CO20180409
  partial_occupation_HNF = bstr.partial_occupation_HNF;
  partial_occupation_sublattice.clear();
  for (size_t i = 0; i < bstr.partial_occupation_sublattice.size(); i++) {
    partial_occupation_sublattice.push_back(bstr.partial_occupation_sublattice[i]);
  }
  // FORCES/POSITIONS ------------------
  qm_calculated = bstr.qm_calculated;
  qm_scale = bstr.qm_scale;
  qm_lattice = bstr.qm_lattice;
  qm_klattice = bstr.qm_klattice;
  qm_f2c = bstr.qm_f2c;
  qm_c2f = bstr.qm_c2f;
  qm_origin = bstr.qm_origin;
  qm_atoms.clear();
  for (size_t i = 0; i < bstr.qm_atoms.size(); i++) {
    qm_atoms.push_back(bstr.qm_atoms[i]);
  }
  qm_forces.clear();
  for (size_t i = 0; i < bstr.qm_forces.size(); i++) {
    qm_forces.push_back(bstr.qm_forces[i]);
  }
  qm_forces_write = bstr.qm_forces_write;
  qm_positions.clear();
  for (size_t i = 0; i < bstr.qm_positions.size(); i++) {
    qm_positions.push_back(bstr.qm_positions[i]);
  }
  qm_positions_write = bstr.qm_positions_write;
  qm_E_cell = bstr.qm_E_cell;
  qm_dE_cell = bstr.qm_dE_cell;
  qm_H_cell = bstr.qm_H_cell;
  qm_PV_cell = bstr.qm_PV_cell;
  qm_mag_cell = bstr.qm_mag_cell;
  qm_P = bstr.qm_P;
  qm_E_atom = bstr.qm_E_atom;
  qm_dE_atom = bstr.qm_dE_atom;
  qm_H_atom = bstr.qm_H_atom;
  qm_PV_atom = bstr.qm_PV_atom;
  qm_mag_atom = bstr.qm_mag_atom;
  // KPOINTS ---------------------------
  kpoints_k1 = bstr.kpoints_k1;
  kpoints_k2 = bstr.kpoints_k2;
  kpoints_k3 = bstr.kpoints_k3;
  kpoints_s1 = bstr.kpoints_s1;
  kpoints_s2 = bstr.kpoints_s2;
  kpoints_s3 = bstr.kpoints_s3;
  kpoints_kmax = bstr.kpoints_kmax;
  kpoints_kppra = bstr.kpoints_kppra;
  kpoints_mode = bstr.kpoints_mode;
  kpoints_kscheme = bstr.kpoints_kscheme;
  // DX+CO START
  dist_nn_min = bstr.dist_nn_min; // CO
  // SYMMETRY TOLERANCE ----------------------------
  sym_eps = bstr.sym_eps; // DX
  sym_eps_calculated = bstr.sym_eps_calculated; // DX
  sym_eps_change_count = bstr.sym_eps_change_count; // DX20180222 - added tolerance count specific to structure
  sym_eps_no_scan = bstr.sym_eps_no_scan; // DX20210331 - added no scan specific to structure
  // DX+CO END
  //  PGROUP ----------------------------
  pgroup.clear();
  for (size_t i = 0; i < bstr.pgroup.size(); i++) {
    pgroup.push_back(bstr.pgroup[i]);
  }
  pgroup_calculated = bstr.pgroup_calculated;
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear();
  for (size_t i = 0; i < bstr.pgroup_xtal.size(); i++) {
    pgroup_xtal.push_back(bstr.pgroup_xtal[i]);
  }
  pgroup_xtal_calculated = bstr.pgroup_xtal_calculated;
  crystal_family = bstr.crystal_family;
  crystal_system = bstr.crystal_system;
  point_group_crystal_class = bstr.point_group_crystal_class;
  point_group_Shoenflies = bstr.point_group_Shoenflies;
  point_group_Hermann_Mauguin = bstr.point_group_Hermann_Mauguin;
  point_group_orbifold = bstr.point_group_orbifold;
  point_group_type = bstr.point_group_type;
  point_group_order = bstr.point_group_order;
  point_group_structure = bstr.point_group_structure;
  // PGROUPK_PATTERSON ---------------------------- //DX20200129
  pgroupk_Patterson.clear();
  for (size_t i = 0; i < bstr.pgroupk_Patterson.size(); i++) {
    pgroupk_Patterson.push_back(bstr.pgroupk_Patterson[i]);
  }
  pgroupk_Patterson_calculated = bstr.pgroupk_Patterson_calculated;
  // PGROUPK ----------------------------
  pgroupk.clear();
  for (size_t i = 0; i < bstr.pgroupk.size(); i++) {
    pgroupk.push_back(bstr.pgroupk[i]);
  }
  pgroupk_calculated = bstr.pgroupk_calculated;
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear(); // DX20171205 - Added pgroupk_xtal
  for (size_t i = 0; i < bstr.pgroupk_xtal.size(); i++) { // DX20171205 - Added pgroupk_xtal
    pgroupk_xtal.push_back(bstr.pgroupk_xtal[i]); // DX20171205 - Added pgroupk_xtal
  }
  pgroupk_xtal_calculated = bstr.pgroupk_xtal_calculated; // DX20171205 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear();
  for (size_t i = 0; i < bstr.fgroup.size(); i++) {
    fgroup.push_back(bstr.fgroup[i]);
  }
  fgroup_calculated = bstr.fgroup_calculated;
  // SGROUP ----------------------------
  sgroup_radius = bstr.sgroup_radius;
  sgroup_radius_dims = bstr.sgroup_radius_dims;
  sgroup.clear();
  for (size_t i = 0; i < bstr.sgroup.size(); i++) {
    sgroup.push_back(bstr.sgroup[i]);
  }
  sgroup_calculated = bstr.sgroup_calculated;
  // SITE POINT GROUP ------------------
  agroup_calculated = bstr.agroup_calculated;
  for (size_t i = 0; i < agroup.size(); i++) {
    agroup[i].clear();
  }
  agroup.clear();
  agroup = std::vector<std::vector<_sym_op>>(bstr.agroup.size());
  for (size_t i = 0; i < bstr.agroup.size(); i++) {
    for (size_t j = 0; j < bstr.agroup[i].size(); j++) {
      agroup[i].push_back(bstr.agroup[i][j]);
    }
  }
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated = bstr.iatoms_calculated;
  for (size_t i = 0; i < iatoms.size(); i++) {
    iatoms[i].clear();
  }
  iatoms.clear();
  for (size_t i = 0; i < bstr.iatoms.size(); i++) {
    iatoms.emplace_back(0);
    for (size_t j = 0; j < bstr.iatoms[i].size(); j++) {
      iatoms.at(i).push_back(bstr.iatoms[i][j]);
    }
  }
  // SPACE GROUP WITH PLATON/FINDSYM/AFLOW -----------
  spacegroup = bstr.spacegroup;
  spacegrouplabel = bstr.spacegrouplabel;
  spacegroupoption = bstr.spacegroupoption;
  spacegroupnumber = bstr.spacegroupnumber;
  spacegroupnumberoption = bstr.spacegroupnumberoption;
  is_spacegroup_platon = bstr.is_spacegroup_platon;
  is_spacegroup_findsym = bstr.is_spacegroup_findsym;
  is_spacegroup_aflow = bstr.is_spacegroup_aflow;
  // SPACE GROUP ITC -----------
  crystal_system_ITC = bstr.crystal_system_ITC; // RHT
  point_group_ITC = bstr.point_group_ITC; // RHT
  bravais_label_ITC = bstr.bravais_label_ITC; // RHT
  lattice_label_ITC = bstr.lattice_label_ITC; // RHT
  space_group_ITC = bstr.space_group_ITC; // RHT
  wyckoff_library_entry_ITC = bstr.wyckoff_library_entry_ITC; // RHT
  wyccar_ITC.clear();
  for (size_t i = 0; i < bstr.wyccar_ITC.size(); i++) {
    wyccar_ITC.push_back(bstr.wyccar_ITC[i]); // RHT
  }
  standard_lattice_ITC = bstr.standard_lattice_ITC; // RHT
  standard_basis_ITC.clear();
  for (size_t i = 0; i < bstr.standard_basis_ITC.size(); i++) {
    standard_basis_ITC.push_back(bstr.standard_basis_ITC[i]); // RHT
  }
  wyckoff_sites_ITC.clear();
  for (size_t i = 0; i < bstr.wyckoff_sites_ITC.size(); i++) {
    wyckoff_sites_ITC.push_back(bstr.wyckoff_sites_ITC[i]);
  }
  // RHT
  wyckoff_symbols_ITC.clear();
  for (size_t i = 0; i < bstr.wyckoff_symbols_ITC.size(); i++) {
    wyckoff_symbols_ITC.push_back(bstr.wyckoff_symbols_ITC[i]); // RHT
  }
  setting_ITC = bstr.setting_ITC; // DX20170830 - SGDATA
  origin_ITC = bstr.origin_ITC; // DX20170830 - SGDATA
  general_position_ITC = bstr.general_position_ITC; // DX20170830 - SGDATA
  // GRID ATOMS ------------------------
  grid_atoms_calculated = bstr.grid_atoms_calculated;
  grid_atoms_dimsL = bstr.grid_atoms_dimsL;
  grid_atoms_dimsH = bstr.grid_atoms_dimsH;
  grid_atoms.clear();
  for (size_t i = 0; i < bstr.grid_atoms.size(); i++) {
    grid_atoms.push_back(bstr.grid_atoms[i]);
  }
  grid_atoms_number = bstr.grid_atoms_number;
  grid_atoms_sc2pcMap.clear();
  for (size_t i = 0; i < bstr.grid_atoms_sc2pcMap.size(); i++) {
    grid_atoms_sc2pcMap.push_back(bstr.grid_atoms_sc2pcMap[i]);
  } // CO20171025
  grid_atoms_pc2scMap.clear();
  for (size_t i = 0; i < bstr.grid_atoms_pc2scMap.size(); i++) {
    grid_atoms_pc2scMap.push_back(bstr.grid_atoms_pc2scMap[i]);
  } // CO20171025
  // LIJK OBEJCTS ----------------------
  lijk_calculated = bstr.lijk_calculated;
  lijk_table.clear();
  lijk_cpos.clear();
  lijk_fpos.clear();
  for (size_t i = 0; i < bstr.lijk_table.size(); i++) {
    lijk_table.push_back(bstr.lijk_table[i]);
    lijk_cpos.push_back(bstr.lijk_cpos.at(i));
    lijk_fpos.push_back(bstr.lijk_fpos.at(i));
  }
  lijk_dims = bstr.lijk_dims;
  // NEIGHBOR ------------------------
  //  for(size_t i=0;i<neighbors_atoms_func_r_vs_nn.size();i++)
  //   neighbors_atoms_func_r_vs_nn.at(i).clear();
  //  for(size_t i=0;i<neighbors_atoms_func_num_vs_nn.size();i++)
  //   neighbors_atoms_func_num_vs_nn.at(i).clear();
  // OUTPUT/ERROR ----------------------
  Niggli_has_failed = bstr.Niggli_has_failed;
  Minkowski_has_failed = bstr.Minkowski_has_failed;
  LatticeReduction_has_failed = bstr.LatticeReduction_has_failed;
  write_lattice_flag = bstr.write_lattice_flag;
  write_klattice_flag = bstr.write_klattice_flag;
  write_inequivalent_flag = bstr.write_inequivalent_flag;
  write_DEBUG_flag = bstr.write_DEBUG_flag;
  error_flag = bstr.error_flag;
  error_string = bstr.error_string;
  // ----------------------------------
}

// ME20200220 - from CO's function in apl::Supercell
void LightCopy(const xstructure& a, xstructure& b) {
  b.clear();
  stringstream POSCAR;
  POSCAR.str("");
  POSCAR << a;
  POSCAR >> b;
  // enable inequivalent flag to work
  for (size_t i = 0; i < b.atoms.size(); i++) {
    b.atoms[i].equivalent = a.atoms[i].equivalent;
    b.atoms[i].is_inequivalent = a.atoms[i].is_inequivalent;
    b.atoms[i].num_equivalents = a.atoms[i].num_equivalents;
  }
  // enable inequivalent flag to work
  b.write_inequivalent_flag = a.write_inequivalent_flag;
  b.info = a.info;
}

// copy
xstructure::xstructure(const xstructure& b) {
  //  free();
  // *this=b;
  copy(b);
}

// destructor
xstructure::~xstructure() {
  free(); // DX20191220 - added free and moved contents below into free
}

// copies xtructures: b=a
const xstructure& xstructure::operator=(const xstructure& b) { // operator=
  if (this != &b) {
    free();
    copy(b);
  }
  return *this;
}

void xstructure::clear() { // DX20191220 - uppercase to lowercase clear
  const xstructure _tmp;
  (*this) = _tmp;
}

void xstructure::clean() { // DX20191220 - uppercase to lowercase clean
  stringstream ss_xstr;
  ss_xstr << (*this);
  (*this).clear(); // DX20191220 - uppercase to lowercase clear
  ss_xstr >> (*this);
  ss_xstr.str("");
}

void xstructure::ClearSpecies() { // CO20180420 - helps with pocc, match with AddAtom()
  num_each_type.clear();
  comp_each_type.clear();
  stoich_each_type.clear();
  species.clear();
  species_pp.clear();
  species_pp_type.clear();
  species_pp_version.clear();
  species_pp_ZVAL.clear();
  species_pp_vLDAU.clear();
  species_volume.clear();
  species_mass.clear();
}

// ME20211004 - from POCC
void xstructure::CleanStructure() {
  neg_scale = false; // NO negative scale... doesn't really matter, scale is one variable
  ReScale(1.0);
  ShiftOriginToAtom(0);
  BringInCell();
  clean(); // DX20191220 - uppercase to lowercase clean
}

void xstructure::initialize(const string& structure_title) {
  // CO20211122 - initialize structure; avoid copying of xstructure
  free(); // DX20191220 - moved contents below into free()
  title = structure_title;
}

void xstructure::initialize(istream& _input, int _iomode) {
  // DX20210129 - initialize structure; avoid copying of xstructure
  free(); // DX20191220 - added free to initialize
  (*this).iomode = _iomode;
  _input >> (*this);
}

void xstructure::initialize(ifstream& _input, int _iomode) {
  // DX20210129 - initialize structure; avoid copying of xstructure
  free(); // DX20191220 - added free to initialize
  (*this).iomode = _iomode;
  _input >> (*this);
}

void xstructure::initialize(const stringstream& __input, int _iomode) {
  // DX20210129 - initialize structure; avoid copying of xstructure
  free(); // DX20191220 - added free to initialize
  (*this).iomode = _iomode;
  stringstream _input(__input.str());
  _input >> (*this);
}

void xstructure::initialize(const string& _input, int _iomode) {
  // CO20211122 - initialize structure; avoid copying of xstructure
  free(); // DX20191220 - added free to initialize
  stringstream strstream;
  aurostd::compressfile2stringstream(_input, strstream); // CO20171025
  (*this).iomode = _iomode;
  (*this).directory = _input; // DX20180526 - location of xstructure
  strstream >> (*this);
}

void xstructure::initialize(const string& url, const string& file, int _iomode) {
  // CO20211122 - initialize structure; avoid copying of xstructure
  free(); // DX20191220 - added free to initialize
  stringstream strstream;
  string url_processed = url;
  if (url.find(":AFLOW") != string::npos) { // HE20220615 safeguard against the direct use of AURLs as suggested by CO
    url_processed = aurostd::StringSubst(url, ":AFLOW", "/AFLOW");
  }
  strstream << aurostd::httpGet(url_processed + "/" + file);
  (*this).iomode = _iomode;
  (*this).directory = url_processed + "/" + file; // DX20180526 - location of xstructure
  strstream >> (*this);
}

// **************************************************************************
// PrintSymbolicMathRepresentation
// **************************************************************************
string xstructure::PrintSymbolicMathRepresentation() {
  xstructure aa((*this));
  const int a_iomode = aa.iomode;
  stringstream oss;
  // VASP
  if (a_iomode == IOAIMS_AUTO || a_iomode == IOVASP_AUTO || a_iomode == IOVASP_POSCAR) { // VASP POSCAR
    stringstream title;
    title << aa.title << " # parameters: " << aa.num_parameters << " # lattice parameters: " << aa.num_lattice_parameters << " # Wyckoff parameters: " << (aa.num_parameters - aa.num_lattice_parameters);
    oss << title.str() << endl;
    oss << "1.0" << endl; // scaling factor
    for (size_t i = 0; i < aa.symbolic_math_lattice.size(); i++) {
      oss << "   " << aurostd::joinWDelimiter(aa.symbolic_math_lattice[i], "  ") << endl;
    }
    if (aa.is_vasp4_poscar_format == true) {} // nothing to do

    if (aa.is_vasp5_poscar_format == true) {
      for (size_t i = 0; i < aa.species.size(); i++) {
        oss << aa.species[i] << " ";
      }
      oss << endl;
    }
    if (aa.partial_occupation_flag == true) {
      for (size_t i = 0; i < aa.num_each_type.size(); i++) {
        oss << aa.num_each_type[i] << "*";
        oss.unsetf(ios_base::floatfield);
        oss << aa.atoms[i].partial_occupation_value << std::fixed << " ";
      }
    } else {
      for (size_t i = 0; i < aa.num_each_type.size(); i++) {
        oss << aa.num_each_type[i] << " ";
      }
    }
    oss << endl;
    if (aa.coord_flag == _COORDS_FRACTIONAL_) {
      oss << "Direct(" << aa.atoms.size() << ") ";
    }
    if (aa.coord_flag == _COORDS_CARTESIAN_) {
      oss << "Cartesian(" << aa.atoms.size() << ") ";
    }
    if (aa.order_parameter_structure == true) {
      oss << "OrderParameter(" << aa.order_parameter_atoms.size() << ") ";
    }
    if (aa.partial_occupation_flag == true) {
      oss << "Partial ";
      // oss.precision(_pocc_precision_);  //CO20170630
      // oss << std::defaultfloat;
      oss.unsetf(ios_base::floatfield);
      oss << "[";
      for (size_t i = 0; i < aa.comp_each_type.size(); i++) {
        oss << char('A' + i) << aa.comp_each_type[i];
      }
      oss << "] ";
      oss << std::fixed;
      // oss.precision(_precision_);       //CO20170630
    } else {
      oss << "[";
      for (size_t i = 0, k = 0; i < aa.num_each_type.size(); k += aa.num_each_type.at(i), i++) {
        oss << char(aa.atoms.at(k).type + 65) << aa.num_each_type.at(i);
      }
      oss << "] ";
    }

    // done
    oss << endl;
    for (size_t i = 0; i < aa.atoms.size(); i++) {
      if (aa.coord_flag == _COORDS_FRACTIONAL_) {
        oss << "   " << aurostd::joinWDelimiter(aa.atoms[i].fpos_equation, "  ");
      } else if (aa.coord_flag == _COORDS_CARTESIAN_) {
        oss << "   " << aurostd::joinWDelimiter(aa.atoms[i].cpos_equation, "  ");
      }
      if (aa.atoms[i].name_is_given == true) {
        oss << " " << aa.atoms[i].name << " ";
        for (uint j = aa.atoms[i].name.length(); j < 5; j++) {
          oss << " ";
        }
      }
      if (aa.partial_occupation_flag == true) {
        oss.unsetf(ios_base::floatfield);
        oss << "pocc=" << aa.atoms[i].partial_occupation_value << "  ";
        oss << std::fixed;
      }
      oss << endl;
    }
  }
  // QE

  // ABINIT

  // AIMS
  if (a_iomode == IOAIMS_AUTO || a_iomode == IOAIMS_GEOM) { // AIMS GEOM
    oss << "# format: symmetry_n_params [n n_lv n_fracpos]" << endl;
    oss << "symmetry_n_params " << aa.num_parameters << " " << aa.num_lattice_parameters << " " << (aa.num_parameters - aa.num_lattice_parameters) << endl;
    // change lattice parameter ratio to separate parameters
    vector<string> parameter_list;
    for (size_t i = 0; i < aa.prototype_parameter_list.size(); i++) {
      if (aa.prototype_parameter_list[i] == "b/a") {
        parameter_list.emplace_back("b");
      } else if (aa.prototype_parameter_list[i] == "c/a") {
        parameter_list.emplace_back("c");
      } else {
        parameter_list.push_back(aa.prototype_parameter_list[i]);
      }
    }
    oss << "symmetry_params " << aurostd::joinWDelimiter(parameter_list, " ") << endl;
    for (uint i = 0; i < 3; i++) {
      oss << "symmetry_lv " << aurostd::joinWDelimiter(aa.symbolic_math_lattice[i], " , ") << endl;
    }
    if (aa.coord_flag == _COORDS_FRACTIONAL_) {
      for (size_t i = 0; i < aa.atoms.size(); i++) {
        oss << "symmetry_frac " << aurostd::joinWDelimiter(aa.atoms[i].fpos_equation, " , ") << endl;
      }
    }
    if (aa.coord_flag == _COORDS_CARTESIAN_) {
      for (size_t i = 0; i < aa.atoms.size(); i++) {
        oss << "symmetry_cart " << aurostd::joinWDelimiter(aa.atoms[i].cpos_equation, " , ") << endl;
      }
    }
  }
  return oss.str();
}

// **************************************************************************
// PrintUNCLE
// **************************************************************************
string xstructure::PrintUNCLE() { // Print in uncle format
  stringstream oss;
  this->iomode = IOVASP_POSCAR;
  oss << "# Structure name:" << endl;
  oss << *this;
  return oss.str();
}

#define cout cout

// **************************************************************************
// xstructure::GetStoich
// **************************************************************************
// Get stoichiometries
bool xstructure::GetStoich() { // CO20171025
  double total_comp = 0.0;
  for (size_t i = 0; i < comp_each_type.size(); i++) {
    total_comp += comp_each_type[i];
  }
  stoich_each_type.clear();
  for (size_t i = 0; i < comp_each_type.size(); i++) {
    stoich_each_type.push_back(comp_each_type[i] / total_comp);
  }
  // CO20210916 - round-off printing errors can be a big challenge here
  // with PARTCAR showing: pocc=1, pocc=0.333, pocc=0.333, pocc=0.333 (LIB4/LIB/CNb_svTa_pvTi_sv:PAW_PBE/AB_cF8_225_a_b.AB:POCC_P0-1xA_P1-0.333xB-0.333xC-0.333xD)
  // we get this:
  // stoichiometry=0.50025013,0.16658329,0.16658329,0.16658329
  // instead of this:
  // stoichiometry=0.5,0.166666667,0.166666667,0.166666667
  // the problem is NOT the partial_occupation_values, but the sum and division for stoich
  // try to fix with aurostd::double2fraction()
  if (partial_occupation_flag) {
    try {
      const bool LDEBUG = (false || XHOST.DEBUG);
      total_comp = 0.0;
      int numerator = 0;
      int denominator = 0;
      double stoich = 0.0;
      deque<double> vstoich;
      for (size_t i = 0; i < comp_each_type.size(); i++) {
        aurostd::double2fraction(comp_each_type[i], numerator, denominator, partial_occupation_stoich_tol);
        // need the right tolerance here
        stoich = (double) numerator / (double) denominator;
        total_comp += stoich;
        vstoich.push_back(stoich);
      }
      for (size_t i = 0; i < comp_each_type.size(); i++) {
        vstoich[i] /= total_comp;
      }
      stoich_each_type = vstoich;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " stoich_each_type=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(stoich_each_type), ",") << endl;
      }
    } catch (aurostd::xerror& re) {
      ;
    } // do nothing
  }
  return true;
}

// **************************************************************************
// xstructure::sortAtomsEquivalent()
// **************************************************************************
// cluster together atoms by equivalent atoms
bool xstructure::sortAtomsEquivalent() {
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (partial_occupation_flag == false) {
    if (!(*this).iatoms_calculated) {
      pflow::PerformFullSymmetry(*this);
    }
    if (!(*this).iatoms_calculated) {
      return false;
    }
  }
  if (LDEBUG) {
    const bool write_inequivalent_flag = (*this).write_inequivalent_flag;
    (*this).write_inequivalent_flag = true;
    cerr << __AFLOW_FUNC__ << " structure before iatoms sorting" << endl;
    cerr << (*this) << endl;
    (*this).write_inequivalent_flag = write_inequivalent_flag;
  }
  bool sort_needed = true;
  if (partial_occupation_flag == false) {
    sort_needed = false;
    int iatom = 0;
    for (size_t i = 0; i < (*this).iatoms.size() && !sort_needed; i++) {
      for (size_t j = 0; j < (*this).iatoms[i].size() && !sort_needed; j++) {
        if ((*this).iatoms[i][j] != iatom) {
          sort_needed = true;
        }
        iatom++;
      }
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " sort " << (sort_needed ? string("") : string("NOT ")) << "needed for this structure!" << endl;
  }
  (*this).MakeBasis(); // make sure the basis is set, sortAtomsEquiv tries not to mess with relative order using basis
  if (!sort_needed) {
    return true;
  }
  deque<_atom> atoms = (*this).atoms;
  std::stable_sort(atoms.begin(), atoms.end(), sortAtomsEquiv); // safe because we do AddAtom() below
  if (LDEBUG) {
    // check order before AddAtom()
    cerr << __AFLOW_FUNC__ << " newly sorted atoms pre-AddAtom()" << endl;
    bool print_RHT = false;
    bool verbose = false;
    bool print_cartesian = false;
    for (size_t i = 0; i < atoms.size(); i++) {
      print_RHT = atoms[i].print_RHT;
      atoms[i].print_RHT = false;
      verbose = atoms[i].verbose;
      atoms[i].verbose = false;
      print_cartesian = atoms[i].print_cartesian;
      atoms[i].print_cartesian = false;
      cerr << atoms[i] << endl;
      atoms[i].print_RHT = print_RHT;
      atoms[i].verbose = verbose;
      atoms[i].print_cartesian = print_cartesian;
    }
  }
  const size_t atoms_size = (*this).atoms.size();
  for (size_t i = atoms_size - 1; i < atoms_size; i--) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " removing atom[" << i << "]" << endl;
    }
    (*this).RemoveAtom(i);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " all atoms removed" << endl;
  }
  for (size_t i = 0; i < atoms.size(); i++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " adding atom[" << i << "]" << endl;
    }
    (*this).AddAtom(atoms[i], false); // CO20230319 - add by type
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " all atoms added back" << endl;
  }
  (*this).ClearSymmetry(); // new structure, clear symmetry
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " newly sorted structure" << endl;
    cerr << (*this) << endl;
  }
  return true;
}

// **************************************************************************
// xstructure::FixLattices
// **************************************************************************
// Fix all the lattices (you can do by hand, but here everything is done)
bool xstructure::FixLattices() {
  klattice = ReciprocalLattice(lattice, scale);
  f2c = trasp(lattice);
  c2f = inverse(trasp(lattice));
  // lattice is the reference... everything else depends
  //  if(iomode==IOVASP_POSCAR || iomode==IOVASP_AUTO) {
  xvector<double> data(6);
  data = Getabc_angles(lattice, DEGREES);
  a = data[1];
  b = data[2];
  c = data[3];
  alpha = data[4];
  beta = data[5];
  gamma = data[6];
  //  }
  //  if(iomode==IOVASP_ABCCAR || iomode==IOVASP_WYCKCAR) {
  //  lattice=GetClat(a,b,c,alpha,beta,gamma);
  // }
  return true;
}

// **************************************************************************
// GetStructure
// **************************************************************************
// get a structure from a directory (POSCAR for VASP
xstructure GetStructure(const int& iomode, ifstream& input) {
  xstructure out;
  if (iomode == IOVASP_POSCAR) { // VASP POSCAR
    out.iomode = IOVASP_POSCAR;
    if (!input) {
      cerr << "EEEEE   File not found" << endl;
      input.clear();
      input.close();
      out.error_flag = true;
      out.error_string = "FILE NOT FOUND";
      return out;
    }
    input >> out;
  }
  return out;
}

// **************************************************************************
// GetStructure
// **************************************************************************
xstructure GetStructure(const int& iomode, const string& Directory) {
  xstructure out;
  if (iomode == IOVASP_POSCAR) { // VASP POSCAR
    out.iomode = IOVASP_POSCAR;
    ifstream input;
    string File;
    File = Directory + "/POSCAR";
    input.open(File.c_str(), std::ios::in);
    out = GetStructure(iomode, input);
    input.clear();
    input.close();
    return out;
  }
  return out;
}

// **************************************************************************
// xstructure::SetCoordinates
// **************************************************************************
// change coordinates type
void xstructure::SetCoordinates(int mode) {
  switch (mode) {
    case _UPDATE_LATTICE_VECTORS_TO_ABCANGLES_: {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "[1] mode=" + aurostd::utype2string(mode), _INPUT_ERROR_);
      break;
    }
    case _UPDATE_LATTICE_ABCANGLES_TO_VECTORS_: {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "[2] mode=" + aurostd::utype2string(mode), _INPUT_ERROR_);
      break;
    }
    case _COORDS_CARTESIAN_: {
      coord_flag = _COORDS_CARTESIAN_;
      strcpy(coord_type, "C");
      break;
    }
    case _COORDS_FRACTIONAL_: {
      coord_flag = _COORDS_FRACTIONAL_;
      strcpy(coord_type, "D");
      break;
    }
    default: {
      cerr << __AFLOW_FUNC__ << " NOTHING TO DO  mode=" << mode << endl;
    }
  }
}

// **************************************************************************
// xstructure::MakeBasis
// **************************************************************************
// This fixes basis and number
void xstructure::MakeBasis() {
  for (size_t iatom = 0; iatom < atoms.size(); iatom++) {
    atoms[iatom].basis = iatom;
    //[CO20200130 - number->basis]atoms[iatom].number=iatom;
  }
}

// **************************************************************************
// xstructure::MakeTypes
// **************************************************************************
// CO20180420
void xstructure::MakeTypes() {
  // need to update TYPES based on num_each_type
  // type is usually used as an index for species
  // if we take a subset of atoms from another structure (POCC), need to reset first iatom to 0
  stringstream message;
  uint sum_atoms = 0;
  for (size_t itype = 0; itype < num_each_type.size(); itype++) {
    sum_atoms += num_each_type[itype];
  }
  if (sum_atoms != atoms.size()) {
    message << "num_each_type does not match atom count (sum_atoms=" << sum_atoms << " vs. atoms.size()=" << atoms.size() << ").";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
  }

  uint iat = 0;
  for (size_t itype = 0; itype < num_each_type.size(); itype++) {
    for (uint j = 0; j < (uint) num_each_type[itype]; j++) {
      atoms[iat++].type = itype;
    }
  }
}

// **************************************************************************
// xstructure::AddAtom() //DX20210202
// **************************************************************************
// This adds a deque<_atom> to the structure.
// More efficient than adding one atom at a time (AddAtom): use more
// efficient for-loop for atoms (upper-triangular) and update species/basis
// info once at the end

void xstructure::AddAtom(const deque<_atom>& _atoms_in, bool add_species, bool check_present) {
  // DX20210129  //CO20230220 - patching names vs. types bug //CO20230319 - adding add_species which assumes the species information of the incoming atoms is correct, otherwise add by type
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    cerr << __AFLOW_FUNC__ << " add_species=" << add_species << endl;
  }

  deque<_atom> atoms_in = _atoms_in; // make a copy so we can fix types

  // CO20230319 = add_species fixes the "source of truth" problem below, we either add by species (and patch type), or we add by type
  // CO20230220 - found a "source of truth" issue: we cannot trust the incoming types of atoms_in,
  // as the type many not match with species of (*this)
  // therefore, we need to rely on the name, and if it's not given, (DO NOT) throw an error, we have to pray that type has been set correctly
  // do NOT modify (*this), as we need to check if the atom sits on top of an existing atom before we really add it to (*this)
  // only modify the type of the incoming atoms vector
  // NOTE: MapAtom() below DEPENDS on having correct type implicitly, so this part of the routine MUST come first
  uint iat = 0;
  if (add_species) { // CO20230319 - fix atom's type
    deque<string> species_new = species;
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ORIG species=" << aurostd::joinWDelimiter(species_new, ",") << endl;
      cerr << __AFLOW_FUNC__ << " NAMES ORIG atoms_in=";
      for (iat = 0; iat < atoms_in.size(); iat++) {
        cerr << atoms_in[iat].name << (iat < atoms_in.size() - 1 ? "," : "");
      }
      cerr << endl;
      cerr << __AFLOW_FUNC__ << " TYPES ORIG atoms_in=";
      for (iat = 0; iat < atoms_in.size(); iat++) {
        cerr << atoms_in[iat].type << (iat < atoms_in.size() - 1 ? "," : "");
      }
      cerr << endl;
    }
    //
    bool FOUND_SPECIES = false;
    uint isp = 0;
    uint species_position = 0;
    for (iat = 0; iat < atoms_in.size(); iat++) {
      _atom& atom = atoms_in[iat];
      if (atom.name.empty()) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "atom.name.empty()", _INPUT_MISSING_);
      }
      FOUND_SPECIES = false;
      for (isp = 0; isp < species_new.size() && FOUND_SPECIES == false; isp++) {
        if (atom.name == species_new[isp]) {
          FOUND_SPECIES = true;
          species_position = isp;
        }
      }
      if (!FOUND_SPECIES) {
        species_new.push_back(atom.name);
        species_position = species_new.size() - 1;
      }
      atom.type = species_position;
    }
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " TYPES NEW atoms_in=";
      for (iat = 0; iat < atoms_in.size(); iat++) {
        cerr << atoms_in[iat].type << (iat < atoms_in.size() - 1 ? "," : "");
      }
      cerr << endl;
      cerr << __AFLOW_FUNC__ << " NEW species=" << aurostd::joinWDelimiter(species_new, ",") << endl;
    }
    //
  }
  for (iat = 0; iat < atoms_in.size(); iat++) {
    _atom& atom = atoms_in[iat];
    atom.CleanName(); // moved up from below
    //[DX20170921 - Need to keep spin info]atom.CleanSpin();
    atom.basis = atoms.size() + iat; // also update basis for sorting below
  }

  const size_t natoms_xstr = atoms.size();
  const deque<_atom>* ptr_atoms = &atoms_in;
  deque<_atom> atoms_unique;

  if (check_present) {
    // use sym_eps if available; if not, use tenth of an Angstrom
    // (since this function adds atoms iteratively, we cannot use minimumDistance,
    // because it would change as we add new atoms) //DX20210202
    double tol = (*this).sym_eps;
    if (tol >= AUROSTD_NAN || tol < _ZERO_TOL_) {
      tol = 0.1;
    } // tenth of Angstrom

    // first check if the input atoms are unique
    // it is more efficient to use a double for-loop (upper-triangular)
    // as opposed to MapAtom(deque<_atom>, _atom); otherwise you end up
    // checking atoms twice //DX20210202
    bool FOUND_POSITION = false;
    uint jat = 0;
    for (iat = 0; iat < atoms_in.size(); iat++) {
      FOUND_POSITION = false;
      for (jat = iat + 1; jat < atoms_in.size() && FOUND_POSITION == false; jat++) {
        if (SYM::MapAtom(atoms_in[iat], atoms_in[jat], true, (*this).lattice, false, tol)) {
          FOUND_POSITION = true;
        }
      }
      if (FOUND_POSITION) {
        continue;
      }
      // now check if any atoms in the xstructure are duplicates with the input atoms
      else if (natoms_xstr != 0) {
        // ME20220120 - changed atoms_unique to atom_in since we are looping over atoms_in
        if (!SYM::MapAtom(atoms, atoms_in[iat], true, (*this).lattice, false, tol)) {
          atoms_unique.push_back(atoms_in[iat]);
        }
      }
      // if no atoms in the xstructure, just add to the unique list
      else {
        atoms_unique.push_back(atoms_in[iat]);
      }
    }
    ptr_atoms = &atoms_unique;
  }

  // update the species: update num/comp each type or add new species
  for (iat = 0; iat < ptr_atoms->size(); iat++) {
    UpdateSpecies(ptr_atoms->at(iat), add_species); // DX20210202 - consolidated code below into function
  }

  // add atoms to xstructure
  if (natoms_xstr == 0) {
    atoms = *ptr_atoms;
  } // if possible, do assignment instead of push_back (faster)
  else {
    for (iat = 0; iat < ptr_atoms->size(); iat++) {
      atoms.push_back(ptr_atoms->at(iat));
    }
  }

  GetStoich(); // CO20170724
  // CO20230220 - BEWARE: this stable sort is NOT meant to sort species
  // it is meant to sort atoms added to the end of the list to the right species-set within already constructed species vector
  // do NOT use sortAtomNames
  //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order]//ME20220612 - Originally had sortAtomsTypes, but other parts of xstructure use sortAtomsNames,
  //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order]//leading to inconsistencies when the input structure is not alphabetic.
  std::stable_sort(atoms.begin(), atoms.end(), sortAtomsTypes);
  //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order] //sortAtomsNames
  MakeBasis(); // need to update NUMBER and BASIS
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " updated xstr=" << endl << (*this) << endl;
  }
}

// **************************************************************************
// xstructure::AddAtom
// **************************************************************************
// This adds an atom to the structure.

void xstructure::AddAtom(const _atom& atom, bool add_species, bool check_present) {
  // CO20230220 - patching names vs. types bug //CO20230319 - adding add_species which assumes the species information of the incoming atoms is correct, otherwise add by type
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    cerr << __AFLOW_FUNC__ << " add_species=" << add_species << endl;
  }

  _atom btom = atom; // DX20210202

  // CO20230319 = add_species fixes the "source of truth" problem below, we either add by species (and patch type), or we add by type
  // CO20230220 - found a "source of truth" issue: we cannot trust the incoming types of atoms_in,
  // as the type many not match with species of (*this)
  // therefore, we need to rely on the name, and if it's not given, (DO NOT) throw an error, we have to pray that type has been set correctly
  // do NOT modify (*this), as we need to check if the atom sits on top of an existing atom before we really add it to (*this)
  // only modify the type of the incoming atoms vector
  // NOTE: MapAtom() below DEPENDS on having correct type implicitly, so this part of the routine MUST come first
  if (add_species) { // CO20230319 - fix atom's type
    if (btom.name.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "atom.name.empty()", _INPUT_MISSING_);
    }
    deque<string> species_new = species;
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ORIG species=" << aurostd::joinWDelimiter(species_new, ",") << endl;
      cerr << __AFLOW_FUNC__ << " TYPES ORIG atom=" << btom.type << endl;
    }
    //
    bool FOUND_SPECIES = false;
    uint species_position = 0;
    uint isp = 0;
    for (isp = 0; isp < species_new.size() && FOUND_SPECIES == false; isp++) {
      if (btom.name == species_new[isp]) {
        FOUND_SPECIES = true;
        species_position = isp;
      }
    }
    if (!FOUND_SPECIES) {
      species_new.push_back(btom.name);
      species_position = species_new.size() - 1;
    }
    btom.type = species_position;
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " TYPES NEW atom=" << btom.type << endl;
      cerr << __AFLOW_FUNC__ << " NEW species=" << aurostd::joinWDelimiter(species_new, ",") << endl;
    }
    //
  }
  btom.CleanName(); // moved up from below
  //[DX20170921 - Need to keep spin info]btom.CleanSpin();
  btom.basis = atoms.size(); // also update basis for sorting below

  if (check_present) { // CO20210116 - AddCorners() should NOT check
    // use sym_eps if available; if not, use tenth of an Angstrom
    // (since this function adds atoms iteratively, we cannot use minimumDistance,
    // because it would change as we add new atoms) //DX20210202
    double tol = (*this).sym_eps;
    if (tol >= AUROSTD_NAN || tol < _ZERO_TOL_) {
      tol = 0.1;
    } // tenth of Angstrom
    if (SYM::MapAtom((*this).atoms, btom, true, (*this).lattice, false, tol)) {
      return;
    }
  }

  // update the species: update num/comp each type or add new species
  UpdateSpecies(btom, add_species); // DX20210202 - consolidated code below into function
  atoms.push_back(btom);
  GetStoich(); // CO20170724
  std::stable_sort(atoms.begin(), atoms.end(), sortAtomsTypes);
  //[CO20230220 - does not work, will rearrange differently than UpdateSpecies() above, so species order does not match atom order] //sortAtomsNames
  MakeBasis(); // need to update NUMBER and BASIS
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " updated xstr=" << endl << (*this) << endl;
  }
}

// **************************************************************************
// xstructure::RemoveAtom
// **************************************************************************
// This removes an atom from the structure.
// CO20170721 - added some safety checks to make sure we weren't deleting an entry
// of a vector/deque that didn't exist (see if size() > itype)
// this is really important because if we build an xstructure on the fly and only occupy
// some of these attributes, the code will seg fault badly, and the error will not
// show up until the xstructure is deconstructed (super confusing)
// make sure to code safely always!
void xstructure::RemoveAtom(const uint& iatom) {
  if (iatom < atoms.size()) {
    const uint itype = atoms.at(iatom).type;
    if (num_each_type.size() > itype) {
      num_each_type.at(itype)--;
    }
    if (comp_each_type.size() > itype) {
      comp_each_type.at(itype) -= atoms.at(iatom).partial_occupation_value;
    }
    if (num_each_type.size() > itype && num_each_type.at(itype) == 0) {
      if (num_each_type.size() > itype) {
        num_each_type.erase(num_each_type.begin() + itype); // erase num_each_type
      }
      if (comp_each_type.size() > itype) {
        comp_each_type.erase(comp_each_type.begin() + itype); // erase comp_each_type
      }
      if (species.size() > itype) {
        species.erase(species.begin() + itype); // erase species
      }
      if (species_pp.size() > itype) {
        species_pp.erase(species_pp.begin() + itype); // erase species_pp
      }
      if (species_pp_type.size() > itype) {
        species_pp_type.erase(species_pp_type.begin() + itype);
      }
      // erase species_pp_type
      if (species_pp_version.size() > itype) {
        species_pp_version.erase(species_pp_version.begin() + itype);
      }
      // erase species_pp_version
      if (species_pp_ZVAL.size() > itype) {
        species_pp_ZVAL.erase(species_pp_ZVAL.begin() + itype);
      }
      // erase species_pp_ZVAL
      if (species_pp_vLDAU.size() > itype) {
        species_pp_vLDAU.erase(species_pp_vLDAU.begin() + itype);
      }
      // erase species_pp_vLDAU
      if (species_volume.size() > itype) {
        species_volume.erase(species_volume.begin() + itype); // erase species_volume
      }
      if (species_mass.size() > itype) {
        species_mass.erase(species_mass.begin() + itype); // erase species_mass
      }
      // CO20170721 - might be better if we wrote function like MakeBasis() for types and did
      // at the end, but it is not unsafe (seg fault) in the way it is written
      for (size_t i = 0; i < atoms.size(); i++) {
        if (i != iatom && atoms[i].type > (int) itype) {
          atoms[i].type--;
        }
      }
    }
    // CO20170721 - this is obsolete with MakeBasis() below!
    //  do atoms
    // for(size_t i=0;i<atoms.size();i++) {
    //   if(i!=iatom && atoms.at(i).number>atoms.at(iatom).number)
    //     atoms.at(i).number--;
    //   if(i!=iatom && atoms.at(i).basis>atoms.at(iatom).basis)
    //     atoms.at(i).basis--;
    // }
    atoms.erase(atoms.begin() + iatom);
    //    // do partial_occupation_flags
    //    for(size_t i=0;i<partial_occupation_flags.size();i++)
    //      if(iatom==partial_occupation_flags.at(i))
    //	partial_occupation_flags.erase(partial_occupation_flags.begin()+i);
    // do order_parameter_atoms
    // CO20170721 - this is okay as it won't seg fault (doesn't delete anything bigger than .size() )
    for (size_t i = 0; i < order_parameter_atoms.size(); i++) {
      if (iatom == order_parameter_atoms[i]) {
        order_parameter_atoms.erase(order_parameter_atoms.begin() + i);
      }
    }
  }
  GetStoich(); // CO20170724
  // done
  MakeBasis(); // need to update NUMBER and BASIS
}

void xstructure::RemoveAtom(vector<uint>& v_atoms_to_remove) { // CO20181226
  const bool LDEBUG = (false || XHOST.DEBUG);
  std::sort(v_atoms_to_remove.begin(), v_atoms_to_remove.end());
  v_atoms_to_remove.erase(std::unique(v_atoms_to_remove.begin(), v_atoms_to_remove.end()), v_atoms_to_remove.end());
  // remove duplicates //CO20181226
  std::sort(v_atoms_to_remove.rbegin(), v_atoms_to_remove.rend());
  // NOTE the r, reverse sort, that way when we remove, it doesn't affect other indices
  for (size_t atom = 0; atom < v_atoms_to_remove.size(); atom++) {
    if (LDEBUG) {
      cerr << "Removing Atom " << v_atoms_to_remove[atom] << endl;
    }
    RemoveAtom(v_atoms_to_remove[atom]);
  }
}

void xstructure::RemoveAtom() { // DX20210129
  // Removes all atoms from an xstructure and clears the related
  // atom/species variable; faster than removing one at a time

  atoms.clear();
  ClearSpecies(); // clears everything species related
  order_parameter_atoms.clear();
}

void xstructure::ReplaceAtoms(const deque<_atom>& new_atoms, bool check_present, bool sort_species) {
  // CO20190520 //DX20210129 - added check_present  //CO20230319 - add bool for sort_species
  // this is the SAFEST/CLEANEST way to replace atoms in an xstructure
  // it takes care of num_each_type, species, etc.
  const bool LDEBUG = (false || XHOST.DEBUG);

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " removing all atoms" << endl;
  }
  RemoveAtom(); // DX20210129 - remove all atoms and clear species variables

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " adding new atoms" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " new_atoms=";
    for (size_t i = 0; i < new_atoms.size(); i++) {
      cerr << new_atoms[i].name << (i < new_atoms.size() - 1 ? "," : " ");
    }
    cerr << endl;
  }
  AddAtom(new_atoms, true, check_present); // adding atoms  //CO20230319 - add by species

  if (sort_species) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PRE-SPECIES-SORT" << endl;
      cerr << __AFLOW_FUNC__ << " xstr=" << endl << (*this) << endl;
      cerr << __AFLOW_FUNC__ << " xstr.species=" << aurostd::joinWDelimiter(species, ",") << endl;
    }

    SpeciesPutAlphabetic(); // DX20210129

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " POST-SPECIES-SORT" << endl;
      cerr << __AFLOW_FUNC__ << " xstr=" << endl << (*this) << endl;
      cerr << __AFLOW_FUNC__ << " xstr.species=" << aurostd::joinWDelimiter(species, ",") << endl;
    }
  }

  // CO20230220 - test that species matches names
  uint itype = 0;
  uint i = 0;
  uint iatom = 0;
  for (itype = 0; itype < num_each_type.size(); itype++) {
    for (i = 0; i < (uint) num_each_type[itype]; i++) {
      if (!atoms[iatom].name.empty() && !species[itype].empty() && atoms[iatom].name != species[itype]) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "atoms and species were incorrectly (re)sorted", _RUNTIME_ERROR_);
      }
      iatom++;
    }
  }
  //
}

// **************************************************************************
// xstructure::RemoveCopies xstructure::RemoveFractionalCopies xstructure::RemoveCartesianCopies
// **************************************************************************
// This removes atoms too close than tol
void xstructure::RemoveCopies(double tol) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  bool flag_isprimitive = false;
  bool flag_remove = false;
  uint iat1 = 0;
  uint iat2 = 0;
  uint irm = 0;
  while (flag_isprimitive == false) {
    if (iat1 == atoms.size()) {
      flag_isprimitive = true;
    } else {
      flag_remove = false;
      for (iat2 = iat1 + 1; iat2 < atoms.size() && flag_remove == false; iat2++) {
        if ((aurostd::modulus(atoms.at(iat1).fpos - atoms.at(iat2).fpos) < tol || aurostd::modulus(atoms.at(iat1).cpos - atoms.at(iat2).cpos) < tol) && flag_remove == false) {
          flag_remove = true;
          irm = iat2;
        }
      }
      if (flag_remove == true) {
        RemoveAtom(irm);
      } else {
        iat1++;
      }
    }
    if (LDEBUG) {
      cout << "DEBUG (RemoveCopies) iat1=" << iat1 << endl;
    }
  }
}

void xstructure::RemoveFractionalCopies(double tol) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  bool flag_isprimitive = false;
  bool flag_remove = false;
  uint iat1 = 0;
  uint iat2 = 0;
  uint irm = 0;
  while (flag_isprimitive == false) {
    if (iat1 == atoms.size()) {
      flag_isprimitive = true;
    } else {
      flag_remove = false;
      for (iat2 = iat1 + 1; iat2 < atoms.size() && flag_remove == false; iat2++) {
        xvector<double> fdiff = atoms.at(iat1).fpos - atoms.at(iat2).fpos; // DX20180503 - account for PBC
        fdiff = SYM::minimizeDistanceFractionalMethod(fdiff); // DX20190613
        if (aurostd::modulus((*this).f2c * fdiff) < tol && flag_remove == false)
        // DX20180503 - account for PBC and perform in Cartesian space
        { // CO20200106 - patching for auto-indenting
          flag_remove = true;
          irm = iat2;
        }
      }
      if (flag_remove == true) {
        RemoveAtom(irm);
      } else {
        iat1++;
      }
    }
    if (LDEBUG) {
      cout << "DEBUG (RemoveFractionalCopies) iat1=" << iat1 << endl;
    }
  }
}

void xstructure::RemoveCartesianCopies(double tol) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  bool flag_isprimitive = false;
  bool flag_remove = false;
  uint iat1 = 0;
  uint iat2 = 0;
  uint irm = 0;
  while (flag_isprimitive == false) {
    if (iat1 == atoms.size()) {
      flag_isprimitive = true;
    } else {
      flag_remove = false;
      for (iat2 = iat1 + 1; iat2 < atoms.size() && flag_remove == false; iat2++) {
        if (aurostd::modulus(atoms.at(iat1).cpos - atoms.at(iat2).cpos) < tol && flag_remove == false) {
          flag_remove = true;
          irm = iat2;
        }
      }
      if (flag_remove == true) {
        RemoveAtom(irm);
      } else {
        iat1++;
      }
    }
    if (LDEBUG) {
      cout << "DEBUG (RemoveCartesianCopies) iat1=" << iat1 << endl;
    }
  }
}

// **************************************************************************
// xstructure::AddCorners
// **************************************************************************
void xstructure::AddCorners() {
  const bool LDEBUG = (false || XHOST.DEBUG);
  xstructure str;
  BringInCell();
  str = *this;
  while (!atoms.empty()) {
    RemoveAtom(0);
  }
  for (size_t iat = 0; iat < str.atoms.size(); iat++) {
    for (double i = 0; i <= 1; i += 0.99) {
      for (double j = 0; j <= 1; j += 0.99) {
        for (double k = 0; k <= 1; k += 0.99) {
          _atom atom = str.atoms[iat];
          atom.fpos[1] += i;
          atom.fpos[2] += j;
          atom.fpos[3] += k;
          atom.cpos = F2C(lattice, atom.fpos);
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " atom.fpos=" << atom.fpos;
          }
          if (atom.fpos[1] <= 1.0 && atom.fpos[2] <= 1.0 && atom.fpos[3] <= 1.0) {
            if (LDEBUG) {
              cerr << " : adding atom";
            }
            AddAtom(atom, false, false);
            // CO20210116 - do NOT check if atom is already there  //CO20230319 - add by type
          }
          if (LDEBUG) {
            cerr << endl;
          }
        }
      }
    }
  }
  title = title + " with_corners";
}

// **************************************************************************
// xstructure::ShiftOriginToAtom
// **************************************************************************
// // Shift the origin to atom(iat)
void xstructure::ShiftOriginToAtom(const int& iat) {
  stringstream message;
  // DX+CO START
  if (iat < 0 || iat >= (int) atoms.size()) {
    message << "iat=" << iat << " out of boundaries (0," << atoms.size() - 1 << ").";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
  }
  xvector<double> frigin(3);
  origin = atoms[iat].cpos;
  frigin = atoms[iat].fpos;
  for (size_t i = 0; i < atoms.size(); i++) {
    atoms[i].fpos = atoms[i].fpos - frigin;
    atoms[i].cpos = atoms[i].cpos - origin;
  }
  // CO TOO SLOW
  // origin=atoms.at(iat).cpos;
  // frigin=atoms.at(iat).fpos;
  // for(size_t i=0;i<atoms.size();i++) {
  //   atoms.at(i).fpos=atoms.at(i).fpos-frigin;
  //   atoms.at(i).cpos=atoms.at(i).cpos-origin;
  // }
  // DX+CO END
}

// **************************************************************************
// SetSDNumbers
// **************************************************************************
xstructure SetSDNumbers(const xstructure& a, const vector<string>& in_sd) {
  // Note that in_sd has one name for each type, not one name for each atom.
  xstructure b(a);
  //  int size=in_sd.size();
  for (size_t cnt = 0; cnt < b.atoms.size(); cnt++) {
    b.atoms[cnt].sd = in_sd[cnt];
    if (in_sd[cnt].size() < 2) {
      cerr << "WARNING:  Atom=" << cnt << " you must specify SD strings 3 characters long (switching to TTT)" << endl;
      b.atoms[cnt].sd = "TTT";
    }
  }
  return b;
}

// **************************************************************************
// SetSDTypes
// **************************************************************************
xstructure SetSDTypes(const xstructure& a, const vector<string>& in_sd) {
  // Note that in_sd has one name for each type, not one name for each atom.
  xstructure b(a);
  int cnt = -1;
  const int size = in_sd.size();
  for (int i = 0; i < size; i++) {
    if (i < (int) b.num_each_type.size()) {
      for (int j = 0; j < (int) b.num_each_type.at(i); j++) {
        cnt++;
        b.atoms.at(cnt).sd = in_sd[i];
        if (in_sd[i].size() < 2) {
          cerr << "WARNING:  Atom=" << cnt << " you must specify SD strings 3 characters long (switching to TTT)" << endl;
          b.atoms.at(cnt).sd = "TTT";
        }
      }
    }
  }
  return b;
}

// **************************************************************************
// GetTypes
// **************************************************************************
vector<int> GetTypes(const xstructure& a) {
  vector<int> out_type;
  for (int i = 0; i < (int) a.atoms.size(); i++) {
    out_type.push_back(a.atoms[i].type);
  }
  return out_type;
}

// **************************************************************************
// GetNames
// **************************************************************************
vector<string> GetNames(const xstructure& a) {
  vector<string> out_name;
  for (int i = 0; i < (int) a.atoms.size(); i++) {
    out_name.push_back(a.atoms[i].name);
  }
  return out_name;
}

// **************************************************************************
// GetCleanNames
// **************************************************************************
vector<string> GetCleanNames(const xstructure& a) {
  vector<string> out_cleanname;
  for (int i = 0; i < (int) a.atoms.size(); i++) {
    out_cleanname.push_back(a.atoms[i].cleanname);
  }
  return out_cleanname;
}

// **************************************************************************
// GetSpins
// **************************************************************************
vector<double> GetSpins(const xstructure& a) {
  vector<double> out_spin;
  for (int i = 0; i < (int) a.atoms.size(); i++) {
    out_spin.push_back(a.atoms[i].spin);
  }
  return out_spin;
}

// ***************************************************************************
// GetElementName
// ***************************************************************************
string GetElementName(string stringin) {
  // need to clean up the _pv stuff of VASP
  for (uint i = 0; i < NUM_ELEMENTS; i++) {
    if (stringin == vatom_symbol.at(i)) {
      return vatom_name.at(i);
    }
  }
  return "NotFound";
}

///@namespace symmetry_lookup
///@brief collection of stmmetry lookup tables
///@author
///@mod{HE,20240531,conversion from old if else structure}
namespace symmetry_lookup {

  // spacegroup number -> spacegroup symbol
  static const std::map<int, std::string> sgn_sg = {
      {  1,               "P1"},
      {  2,              "P-1"},
      {  3,               "P2"},
      {  4,           "P2_{1}"},
      {  5,               "C2"},
      {  6,               "Pm"},
      {  7,               "Pc"},
      {  8,               "Cm"},
      {  9,               "Cc"},
      { 10,             "P2/m"},
      { 11,         "P2_{1}/m"},
      { 12,             "C2/m"},
      { 13,             "P2/c"},
      { 14,         "P2_{1}/c"},
      { 15,             "C2/c"},
      { 16,             "P222"},
      { 17,         "P222_{1}"},
      { 18,     "P2_{1}2_{1}2"},
      { 19, "P2_{1}2_{1}2_{1}"},
      { 20,         "C222_{1}"},
      { 21,             "C222"},
      { 22,             "F222"},
      { 23,             "I222"},
      { 24, "I2_{1}2_{1}2_{1}"},
      { 25,             "Pmm2"},
      { 26,         "Pmc2_{1}"},
      { 27,             "Pcc2"},
      { 28,             "Pma2"},
      { 29,         "Pca2_{1}"},
      { 30,             "Pnc2"},
      { 31,         "Pmn2_{1}"},
      { 32,             "Pba2"},
      { 33,         "Pna2_{1}"},
      { 34,             "Pnn2"},
      { 35,             "Cmm2"},
      { 36,         "Cmc2_{1}"},
      { 37,             "Ccc2"},
      { 38,             "Amm2"},
      { 39,             "Aem2"},
      { 40,             "Ama2"},
      { 41,             "Aea2"},
      { 42,             "Fmm2"},
      { 43,             "Fdd2"},
      { 44,             "Imm2"},
      { 45,             "Iba2"},
      { 46,             "Ima2"},
      { 47,             "Pmmm"},
      { 48,             "Pnnn"},
      { 49,             "Pccm"},
      { 50,             "Pban"},
      { 51,             "Pmma"},
      { 52,             "Pnna"},
      { 53,             "Pmna"},
      { 54,             "Pcca"},
      { 55,             "Pbam"},
      { 56,             "Pccn"},
      { 57,             "Pbcm"},
      { 58,             "Pnnm"},
      { 59,             "Pmmn"},
      { 60,             "Pbcn"},
      { 61,             "Pbca"},
      { 62,             "Pnma"},
      { 63,             "Cmcm"},
      { 64,             "Cmce"},
      { 65,             "Cmmm"},
      { 66,             "Cccm"},
      { 67,             "Cmme"},
      { 68,             "Ccce"},
      { 69,             "Fmmm"},
      { 70,             "Fddd"},
      { 71,             "Immm"},
      { 72,             "Ibam"},
      { 73,             "Ibca"},
      { 74,             "Imma"},
      { 75,               "P4"},
      { 76,           "P4_{1}"},
      { 77,           "P4_{2}"},
      { 78,           "P4_{3}"},
      { 79,               "I4"},
      { 80,           "I4_{1}"},
      { 81,              "P-4"},
      { 82,              "I-4"},
      { 83,             "P4/m"},
      { 84,         "P4_{2}/m"},
      { 85,             "P4/n"},
      { 86,         "P4_{2}/n"},
      { 87,             "I4/m"},
      { 88,         "I4_{1}/a"},
      { 89,             "P422"},
      { 90,         "P42_{1}2"},
      { 91,         "P4_{1}22"},
      { 92,     "P4_{1}2_{1}2"},
      { 93,         "P4_{2}22"},
      { 94,     "P4_{2}2_{1}2"},
      { 95,         "P4_{3}22"},
      { 96,     "P4_{3}2_{1}2"},
      { 97,             "I422"},
      { 98,         "I4_{1}22"},
      { 99,             "P4mm"},
      {100,             "P4bm"},
      {101,         "P4_{2}cm"},
      {102,         "P4_{2}nm"},
      {103,             "P4cc"},
      {104,             "P4nc"},
      {105,         "P4_{2}mc"},
      {106,         "P4_{2}bc"},
      {107,             "I4mm"},
      {108,             "I4cm"},
      {109,         "I4_{1}md"},
      {110,         "I4_{1}cd"},
      {111,            "P-42m"},
      {112,            "P-42c"},
      {113,        "P-42_{1}m"},
      {114,        "P-42_{1}c"},
      {115,            "P-4m2"},
      {116,            "P-4c2"},
      {117,            "P-4b2"},
      {118,            "P-4n2"},
      {119,            "I-4m2"},
      {120,            "I-4c2"},
      {121,            "I-42m"},
      {122,            "I-42d"},
      {123,           "P4/mmm"},
      {124,           "P4/mcc"},
      {125,           "P4/nbm"},
      {126,           "P4/nnc"},
      {127,           "P4/mbm"},
      {128,           "P4/mnc"},
      {129,           "P4/nmm"},
      {130,           "P4/ncc"},
      {131,       "P4_{2}/mmc"},
      {132,       "P4_{2}/mcm"},
      {133,       "P4_{2}/nbc"},
      {134,       "P4_{2}/nnm"},
      {135,       "P4_{2}/mbc"},
      {136,       "P4_{2}/mnm"},
      {137,       "P4_{2}/nmc"},
      {138,       "P4_{2}/ncm"},
      {139,           "I4/mmm"},
      {140,           "I4/mcm"},
      {141,       "I4_{1}/amd"},
      {142,       "I4_{1}/acd"},
      {143,               "P3"},
      {144,           "P3_{1}"},
      {145,           "P3_{2}"},
      {146,               "R3"},
      {147,              "P-3"},
      {148,              "R-3"},
      {149,             "P312"},
      {150,             "P321"},
      {151,         "P3_{1}12"},
      {152,         "P3_{1}21"},
      {153,         "P3_{2}12"},
      {154,         "P3_{2}21"},
      {155,              "R32"},
      {156,             "P3m1"},
      {157,             "P31m"},
      {158,             "P3c1"},
      {159,             "P31c"},
      {160,              "R3m"},
      {161,              "R3c"},
      {162,            "P-31m"},
      {163,            "P-31c"},
      {164,            "P-3m1"},
      {165,            "P-3c1"},
      {166,             "R-3m"},
      {167,             "R-3c"},
      {168,               "P6"},
      {169,           "P6_{1}"},
      {170,           "P6_{5}"},
      {171,           "P6_{2}"},
      {172,           "P6_{4}"},
      {173,           "P6_{3}"},
      {174,              "P-6"},
      {175,             "P6/m"},
      {176,         "P6_{3}/m"},
      {177,             "P622"},
      {178,         "P6_{1}22"},
      {179,         "P6_{5}22"},
      {180,         "P6_{2}22"},
      {181,         "P6_{4}22"},
      {182,         "P6_{3}22"},
      {183,             "P6mm"},
      {184,             "P6cc"},
      {185,         "P6_{3}cm"},
      {186,         "P6_{3}mc"},
      {187,            "P-6m2"},
      {188,            "P-6c2"},
      {189,            "P-62m"},
      {190,            "P-62c"},
      {191,           "P6/mmm"},
      {192,           "P6/mcc"},
      {193,       "P6_{3}/mcm"},
      {194,       "P6_{3}/mmc"},
      {195,              "P23"},
      {196,              "F23"},
      {197,              "I23"},
      {198,          "P2_{1}3"},
      {199,          "I2_{1}3"},
      {200,             "Pm-3"},
      {201,             "Pn-3"},
      {202,             "Fm-3"},
      {203,             "Fd-3"},
      {204,             "Im-3"},
      {205,             "Pa-3"},
      {206,             "Ia-3"},
      {207,             "P432"},
      {208,         "P4_{2}32"},
      {209,             "F432"},
      {210,         "F4_{1}32"},
      {211,             "I432"},
      {212,         "P4_{3}32"},
      {213,         "P4_{1}32"},
      {214,         "I4_{1}32"},
      {215,            "P-43m"},
      {216,            "F-43m"},
      {217,            "I-43m"},
      {218,            "P-43n"},
      {219,            "F-43c"},
      {220,            "I-43d"},
      {221,            "Pm-3m"},
      {222,            "Pn-3n"},
      {223,            "Pm-3n"},
      {224,            "Pn-3m"},
      {225,            "Fm-3m"},
      {226,            "Fm-3c"},
      {227,            "Fd-3m"},
      {228,            "Fd-3c"},
      {229,            "Im-3m"},
      {230,            "Ia-3d"}
  }; ///< spacegroup number -> spacegroup symbol

  // spacegroup symbol -> spacegroup number
  static const std::map<std::string, int> sg_sgn = {
      {              "P1",   1},
      {             "P-1",   2},
      {              "P2",   3},
      {          "P2_{1}",   4},
      {            "P2_1",   4},
      {             "P21",   4},
      {              "C2",   5},
      {              "Pm",   6},
      {              "Pc",   7},
      {              "Cm",   8},
      {              "Cc",   9},
      {            "P2/m",  10},
      {        "P2_{1}/m",  11},
      {          "P2_1/m",  11},
      {           "P21/m",  11},
      {            "C2/m",  12},
      {            "P2/c",  13},
      {        "P2_{1}/c",  14},
      {          "P2_1/c",  14},
      {           "P21/c",  14},
      {            "C2/c",  15},
      {            "P222",  16},
      {        "P222_{1}",  17},
      {          "P222_1",  17},
      {           "P2221",  17},
      {    "P2_{1}2_{1}2",  18},
      {        "P2_12_12",  18},
      {          "P21212",  18},
      {"P2_{1}2_{1}2_{1}",  19},
      {      "P2_12_12_1",  19},
      {         "P212121",  19},
      {        "C222_{1}",  20},
      {          "C222_1",  20},
      {           "C2221",  20},
      {            "C222",  21},
      {            "F222",  22},
      {            "I222",  23},
      {"I2_{1}2_{1}2_{1}",  24},
      {      "I2_12_12_1",  24},
      {         "I212121",  24},
      {            "Pmm2",  25},
      {        "Pmc2_{1}",  26},
      {          "Pmc2_1",  26},
      {           "Pmc21",  26},
      {            "Pcc2",  27},
      {            "Pma2",  28},
      {        "Pca2_{1}",  29},
      {          "Pca2_1",  29},
      {           "Pca21",  29},
      {            "Pnc2",  30},
      {        "Pmn2_{1}",  31},
      {          "Pmn2_1",  31},
      {           "Pmn21",  31},
      {            "Pba2",  32},
      {        "Pna2_{1}",  33},
      {          "Pna2_1",  33},
      {           "Pna21",  33},
      {            "Pnn2",  34},
      {            "Cmm2",  35},
      {        "Cmc2_{1}",  36},
      {          "Cmc2_1",  36},
      {           "Cmc21",  36},
      {            "Ccc2",  37},
      {            "Amm2",  38},
      {            "Aem2",  39},
      {            "Ama2",  40},
      {            "Aea2",  41},
      {            "Fmm2",  42},
      {            "Fdd2",  43},
      {            "Imm2",  44},
      {            "Iba2",  45},
      {            "Ima2",  46},
      {            "Pmmm",  47},
      {            "Pnnn",  48},
      {            "Pccm",  49},
      {            "Pban",  50},
      {            "Pmma",  51},
      {            "Pnna",  52},
      {            "Pmna",  53},
      {            "Pcca",  54},
      {            "Pbam",  55},
      {            "Pccn",  56},
      {            "Pbcm",  57},
      {            "Pnnm",  58},
      {            "Pmmn",  59},
      {            "Pbcn",  60},
      {            "Pbca",  61},
      {            "Pnma",  62},
      {            "Cmcm",  63},
      {            "Cmce",  64},
      {            "Cmmm",  65},
      {            "Cccm",  66},
      {            "Cmme",  67},
      {            "Ccce",  68},
      {            "Fmmm",  69},
      {            "Fddd",  70},
      {            "Immm",  71},
      {            "Ibam",  72},
      {            "Ibca",  73},
      {            "Imma",  74},
      {              "P4",  75},
      {          "P4_{1}",  76},
      {            "P4_1",  76},
      {             "P41",  76},
      {          "P4_{2}",  77},
      {            "P4_2",  77},
      {             "P42",  77},
      {          "P4_{3}",  78},
      {            "P4_3",  78},
      {             "P43",  78},
      {              "I4",  79},
      {          "I4_{1}",  80},
      {            "I4_1",  80},
      {             "I41",  80},
      {             "P-4",  81},
      {             "I-4",  82},
      {            "P4/m",  83},
      {        "P4_{2}/m",  84},
      {          "P4_2/m",  84},
      {           "P42/m",  84},
      {            "P4/n",  85},
      {        "P4_{2}/n",  86},
      {          "P4_2/n",  86},
      {           "P42/n",  86},
      {            "I4/m",  87},
      {        "I4_{1}/a",  88},
      {          "I4_1/a",  88},
      {           "I41/a",  88},
      {            "P422",  89},
      {        "P42_{1}2",  90},
      {          "P42_12",  90},
      {           "P4212",  90},
      {        "P4_{1}22",  91},
      {          "P4_122",  91},
      {           "P4122",  91},
      {    "P4_{1}2_{1}2",  92},
      {        "P4_12_12",  92},
      {          "P41212",  92},
      {        "P4_{2}22",  93},
      {          "P4_222",  93},
      {           "P4222",  93},
      {    "P4_{2}2_{1}2",  94},
      {        "P4_22_12",  94},
      {          "P42212",  94},
      {        "P4_{3}22",  95},
      {          "P4_322",  95},
      {           "P4322",  95},
      {    "P4_{3}2_{1}2",  96},
      {        "P4_32_12",  96},
      {          "P43212",  96},
      {            "I422",  97},
      {        "I4_{1}22",  98},
      {          "I4_122",  98},
      {           "I4122",  98},
      {            "P4mm",  99},
      {            "P4bm", 100},
      {        "P4_{2}cm", 101},
      {          "P4_2cm", 101},
      {           "P42cm", 101},
      {        "P4_{2}nm", 102},
      {          "P4_2nm", 102},
      {           "P42nm", 102},
      {            "P4cc", 103},
      {            "P4nc", 104},
      {        "P4_{2}mc", 105},
      {          "P4_2mc", 105},
      {           "P42mc", 105},
      {        "P4_{2}bc", 106},
      {          "P4_2bc", 106},
      {           "P42bc", 106},
      {            "I4mm", 107},
      {            "I4cm", 108},
      {        "I4_{1}md", 109},
      {          "I4_1md", 109},
      {           "I41md", 109},
      {        "I4_{1}cd", 110},
      {          "I4_1cd", 110},
      {           "I41cd", 110},
      {           "P-42m", 111},
      {           "P-42c", 112},
      {       "P-42_{1}m", 113},
      {         "P-42_1m", 113},
      {          "P-421m", 113},
      {       "P-42_{1}c", 114},
      {         "P-42_1c", 114},
      {          "P-421c", 114},
      {           "P-4m2", 115},
      {           "P-4c2", 116},
      {           "P-4b2", 117},
      {           "P-4n2", 118},
      {           "I-4m2", 119},
      {           "I-4c2", 120},
      {           "I-42m", 121},
      {           "I-42d", 122},
      {          "P4/mmm", 123},
      {          "P4/mcc", 124},
      {          "P4/nbm", 125},
      {          "P4/nnc", 126},
      {          "P4/mbm", 127},
      {          "P4/mnc", 128},
      {          "P4/nmm", 129},
      {          "P4/ncc", 130},
      {      "P4_{2}/mmc", 131},
      {        "P4_2/mmc", 131},
      {         "P42/mmc", 131},
      {      "P4_{2}/mcm", 132},
      {        "P4_2/mcm", 132},
      {         "P42/mcm", 132},
      {      "P4_{2}/nbc", 133},
      {        "P4_2/nbc", 133},
      {         "P42/nbc", 133},
      {      "P4_{2}/nnm", 134},
      {        "P4_2/nnm", 134},
      {         "P42/nnm", 134},
      {      "P4_{2}/mbc", 135},
      {        "P4_2/mbc", 135},
      {         "P42/mbc", 135},
      {      "P4_{2}/mnm", 136},
      {        "P4_2/mnm", 136},
      {         "P42/mnm", 136},
      {      "P4_{2}/nmc", 137},
      {        "P4_2/nmc", 137},
      {         "P42/nmc", 137},
      {      "P4_{2}/ncm", 138},
      {        "P4_2/ncm", 138},
      {         "P42/ncm", 138},
      {          "I4/mmm", 139},
      {          "I4/mcm", 140},
      {      "I4_{1}/amd", 141},
      {        "I4_1/amd", 141},
      {         "I41/amd", 141},
      {      "I4_{1}/acd", 142},
      {        "I4_1/acd", 142},
      {         "I41/acd", 142},
      {              "P3", 143},
      {          "P3_{1}", 144},
      {            "P3_1", 144},
      {             "P31", 144},
      {          "P3_{2}", 145},
      {            "P3_2", 145},
      {             "P32", 145},
      {              "R3", 146},
      {             "P-3", 147},
      {             "R-3", 148},
      {            "P312", 149},
      {            "P321", 150},
      {        "P3_{1}12", 151},
      {          "P3_112", 151},
      {           "P3112", 151},
      {        "P3_{1}21", 152},
      {          "P3_121", 152},
      {           "P3121", 152},
      {        "P3_{2}12", 153},
      {          "P3_212", 153},
      {           "P3212", 153},
      {        "P3_{2}21", 154},
      {          "P3_221", 154},
      {           "P3221", 154},
      {             "R32", 155},
      {            "P3m1", 156},
      {            "P31m", 157},
      {            "P3c1", 158},
      {            "P31c", 159},
      {             "R3m", 160},
      {             "R3c", 161},
      {           "P-31m", 162},
      {           "P-31c", 163},
      {           "P-3m1", 164},
      {           "P-3c1", 165},
      {            "R-3m", 166},
      {            "R-3c", 167},
      {              "P6", 168},
      {          "P6_{1}", 169},
      {            "P6_1", 169},
      {             "P61", 169},
      {          "P6_{5}", 170},
      {            "P6_5", 170},
      {             "P65", 170},
      {          "P6_{2}", 171},
      {            "P6_2", 171},
      {             "P62", 171},
      {          "P6_{4}", 172},
      {            "P6_4", 172},
      {             "P64", 172},
      {          "P6_{3}", 173},
      {            "P6_3", 173},
      {             "P63", 173},
      {             "P-6", 174},
      {            "P6/m", 175},
      {        "P6_{3}/m", 176},
      {          "P6_3/m", 176},
      {           "P63/m", 176},
      {            "P622", 177},
      {        "P6_{1}22", 178},
      {          "P6_122", 178},
      {           "P6122", 178},
      {        "P6_{5}22", 179},
      {          "P6_522", 179},
      {           "P6522", 179},
      {        "P6_{2}22", 180},
      {          "P6_222", 180},
      {           "P6222", 180},
      {        "P6_{4}22", 181},
      {          "P6_422", 181},
      {           "P6422", 181},
      {        "P6_{3}22", 182},
      {          "P6_322", 182},
      {           "P6322", 182},
      {            "P6mm", 183},
      {            "P6cc", 184},
      {        "P6_{3}cm", 185},
      {          "P6_3cm", 185},
      {           "P63cm", 185},
      {        "P6_{3}mc", 186},
      {          "P6_3mc", 186},
      {           "P63mc", 186},
      {           "P-6m2", 187},
      {           "P-6c2", 188},
      {           "P-62m", 189},
      {           "P-62c", 190},
      {          "P6/mmm", 191},
      {          "P6/mcc", 192},
      {      "P6_{3}/mcm", 193},
      {        "P6_3/mcm", 193},
      {         "P63/mcm", 193},
      {      "P6_{3}/mmc", 194},
      {        "P6_3/mmc", 194},
      {         "P63/mmc", 194},
      {             "P23", 195},
      {             "F23", 196},
      {             "I23", 197},
      {         "P2_{1}3", 198},
      {           "P2_13", 198},
      {            "P213", 198},
      {         "I2_{1}3", 199},
      {           "I2_13", 199},
      {            "I213", 199},
      {            "Pm-3", 200},
      {            "Pn-3", 201},
      {            "Fm-3", 202},
      {            "Fd-3", 203},
      {            "Im-3", 204},
      {            "Pa-3", 205},
      {            "Ia-3", 206},
      {            "P432", 207},
      {        "P4_{2}32", 208},
      {          "P4_232", 208},
      {           "P4232", 208},
      {            "F432", 209},
      {        "F4_{1}32", 210},
      {          "F4_132", 210},
      {           "F4132", 210},
      {            "I432", 211},
      {        "P4_{3}32", 212},
      {          "P4_332", 212},
      {           "P4332", 212},
      {        "P4_{1}32", 213},
      {          "P4_132", 213},
      {           "P4132", 213},
      {        "I4_{1}32", 214},
      {          "I4_132", 214},
      {           "I4132", 214},
      {           "P-43m", 215},
      {           "F-43m", 216},
      {           "I-43m", 217},
      {           "P-43n", 218},
      {           "F-43c", 219},
      {           "I-43d", 220},
      {           "Pm-3m", 221},
      {           "Pn-3n", 222},
      {           "Pm-3n", 223},
      {           "Pn-3m", 224},
      {           "Fm-3m", 225},
      {           "Fm-3c", 226},
      {           "Fd-3m", 227},
      {           "Fd-3c", 228},
      {           "Im-3m", 229},
      {           "Ia-3d", 230}
  }; ///< spacegroup symbol -> spacegroup number

  // spacegroup number -> Schoenflies notation
  static const std::map<int, std::string> sgn_sch = {
      {  1,   "C_{1}^{1}"},
      {  2,   "C_{i}^{1}"},
      {  3,   "C_{2}^{1}"},
      {  4,   "C_{2}^{2}"},
      {  5,   "C_{2}^{3}"},
      {  6,   "C_{s}^{1}"},
      {  7,   "C_{s}^{2}"},
      {  8,   "C_{s}^{3}"},
      {  9,   "C_{s}^{4}"},
      { 10,  "C_{2h}^{1}"},
      { 11,  "C_{2h}^{2}"},
      { 12,  "C_{2h}^{3}"},
      { 13,  "C_{2h}^{4}"},
      { 14,  "C_{2h}^{5}"},
      { 15,  "C_{2h}^{6}"},
      { 16,   "D_{2}^{1}"},
      { 17,   "D_{2}^{2}"},
      { 18,   "D_{2}^{3}"},
      { 19,   "D_{2}^{4}"},
      { 20,   "D_{2}^{5}"},
      { 21,   "D_{2}^{6}"},
      { 22,   "D_{2}^{7}"},
      { 23,   "D_{2}^{8}"},
      { 24,   "D_{2}^{9}"},
      { 25,  "C_{2v}^{1}"},
      { 26,  "C_{2v}^{2}"},
      { 27,  "C_{2v}^{3}"},
      { 28,  "C_{2v}^{4}"},
      { 29,  "C_{2v}^{5}"},
      { 30,  "C_{2v}^{6}"},
      { 31,  "C_{2v}^{7}"},
      { 32,  "C_{2v}^{8}"},
      { 33,  "C_{2v}^{9}"},
      { 34, "C_{2v}^{10}"},
      { 35, "C_{2v}^{11}"},
      { 36, "C_{2v}^{12}"},
      { 37, "C_{2v}^{13}"},
      { 38, "C_{2v}^{14}"},
      { 39, "C_{2v}^{15}"},
      { 40, "C_{2v}^{16}"},
      { 41, "C_{2v}^{17}"},
      { 42, "C_{2v}^{18}"},
      { 43, "C_{2v}^{19}"},
      { 44, "C_{2v}^{20}"},
      { 45, "C_{2v}^{21}"},
      { 46, "C_{2v}^{22}"},
      { 47,  "D_{2h}^{1}"},
      { 48,  "D_{2h}^{2}"},
      { 49,  "D_{2h}^{3}"},
      { 50,  "D_{2h}^{4}"},
      { 51,  "D_{2h}^{5}"},
      { 52,  "D_{2h}^{6}"},
      { 53,  "D_{2h}^{7}"},
      { 54,  "D_{2h}^{8}"},
      { 55,  "D_{2h}^{9}"},
      { 56, "D_{2h}^{10}"},
      { 57, "D_{2h}^{11}"},
      { 58, "D_{2h}^{12}"},
      { 59, "D_{2h}^{13}"},
      { 60, "D_{2h}^{14}"},
      { 61, "D_{2h}^{15}"},
      { 62, "D_{2h}^{16}"},
      { 63, "D_{2h}^{17}"},
      { 64, "D_{2h}^{18}"},
      { 65, "D_{2h}^{19}"},
      { 66, "D_{2h}^{20}"},
      { 67, "D_{2h}^{21}"},
      { 68, "D_{2h}^{22}"},
      { 69, "D_{2h}^{23}"},
      { 70, "D_{2h}^{24}"},
      { 71, "D_{2h}^{25}"},
      { 72, "D_{2h}^{26}"},
      { 73, "D_{2h}^{27}"},
      { 74, "D_{2h}^{28}"},
      { 75,   "C_{4}^{1}"},
      { 76,   "C_{4}^{2}"},
      { 77,   "C_{4}^{3}"},
      { 78,   "C_{4}^{4}"},
      { 79,   "C_{4}^{5}"},
      { 80,   "C_{4}^{6}"},
      { 81,   "S_{4}^{1}"},
      { 82,   "S_{4}^{2}"},
      { 83,  "C_{4h}^{1}"},
      { 84,  "C_{4h}^{2}"},
      { 85,  "C_{4h}^{3}"},
      { 86,  "C_{4h}^{4}"},
      { 87,  "C_{4h}^{5}"},
      { 88,  "C_{4h}^{6}"},
      { 89,   "D_{4}^{1}"},
      { 90,   "D_{4}^{2}"},
      { 91,   "D_{4}^{3}"},
      { 92,   "D_{4}^{4}"},
      { 93,   "D_{4}^{5}"},
      { 94,   "D_{4}^{6}"},
      { 95,   "D_{4}^{7}"},
      { 96,   "D_{4}^{8}"},
      { 97,   "D_{4}^{9}"},
      { 98,  "D_{4}^{10}"},
      { 99,  "C_{4v}^{1}"},
      {100,  "C_{4v}^{2}"},
      {101,  "C_{4v}^{3}"},
      {102,  "C_{4v}^{4}"},
      {103,  "C_{4v}^{5}"},
      {104,  "C_{4v}^{6}"},
      {105,  "C_{4v}^{7}"},
      {106,  "C_{4v}^{8}"},
      {107,  "C_{4v}^{9}"},
      {108, "C_{4v}^{10}"},
      {109, "C_{4v}^{11}"},
      {110, "C_{4v}^{12}"},
      {111,  "D_{2d}^{1}"},
      {112,  "D_{2d}^{2}"},
      {113,  "D_{2d}^{3}"},
      {114,  "D_{2d}^{4}"},
      {115,  "D_{2d}^{5}"},
      {116,  "D_{2d}^{6}"},
      {117,  "D_{2d}^{7}"},
      {118,  "D_{2d}^{8}"},
      {119,  "D_{2d}^{9}"},
      {120, "D_{2d}^{10}"},
      {121, "D_{2d}^{11}"},
      {122, "D_{2d}^{12}"},
      {123,  "D_{4h}^{1}"},
      {124,  "D_{4h}^{2}"},
      {125,  "D_{4h}^{3}"},
      {126,  "D_{4h}^{4}"},
      {127,  "D_{4h}^{5}"},
      {128,  "D_{4h}^{6}"},
      {129,  "D_{4h}^{7}"},
      {130,  "D_{4h}^{8}"},
      {131,  "D_{4h}^{9}"},
      {132, "D_{4h}^{10}"},
      {133, "D_{4h}^{11}"},
      {134, "D_{4h}^{12}"},
      {135, "D_{4h}^{13}"},
      {136, "D_{4h}^{14}"},
      {137, "D_{4h}^{15}"},
      {138, "D_{4h}^{16}"},
      {139, "D_{4h}^{17}"},
      {140, "D_{4h}^{18}"},
      {141, "D_{4h}^{19}"},
      {142, "D_{4h}^{20}"},
      {143,   "C_{3}^{1}"},
      {144,   "C_{3}^{2}"},
      {145,   "C_{3}^{3}"},
      {146,   "C_{3}^{4}"},
      {147,  "C_{3i}^{1}"},
      {148,  "C_{3i}^{2}"},
      {149,   "D_{3}^{1}"},
      {150,   "D_{3}^{2}"},
      {151,   "D_{3}^{3}"},
      {152,   "D_{3}^{4}"},
      {153,   "D_{3}^{5}"},
      {154,   "D_{3}^{6}"},
      {155,   "D_{3}^{7}"},
      {156,  "C_{3v}^{1}"},
      {157,  "C_{3v}^{2}"},
      {158,  "C_{3v}^{3}"},
      {159,  "C_{3v}^{4}"},
      {160,  "C_{3v}^{5}"},
      {161,  "C_{3v}^{6}"},
      {162,  "D_{3d}^{1}"},
      {163,  "D_{3d}^{2}"},
      {164,  "D_{3d}^{3}"},
      {165,  "D_{3d}^{4}"},
      {166,  "D_{3d}^{5}"},
      {167,  "D_{3d}^{6}"},
      {168,   "C_{6}^{1}"},
      {169,   "C_{6}^{2}"},
      {170,   "C_{6}^{3}"},
      {171,   "C_{6}^{4}"},
      {172,   "C_{6}^{5}"},
      {173,   "C_{6}^{6}"},
      {174,  "C_{3h}^{1}"},
      {175,  "C_{6h}^{1}"},
      {176,  "C_{6h}^{2}"},
      {177,   "D_{6}^{1}"},
      {178,   "D_{6}^{2}"},
      {179,   "D_{6}^{3}"},
      {180,   "D_{6}^{4}"},
      {181,   "D_{6}^{5}"},
      {182,   "D_{6}^{6}"},
      {183,  "C_{6v}^{1}"},
      {184,  "C_{6v}^{2}"},
      {185,  "C_{6v}^{3}"},
      {186,  "C_{6v}^{4}"},
      {187,  "D_{3h}^{1}"},
      {188,  "D_{3h}^{2}"},
      {189,  "D_{3h}^{3}"},
      {190,  "D_{3h}^{4}"},
      {191,  "D_{6h}^{1}"},
      {192,  "D_{6h}^{2}"},
      {193,  "D_{6h}^{3}"},
      {194,  "D_{6h}^{4}"},
      {195,       "T^{1}"},
      {196,       "T^{2}"},
      {197,       "T^{3}"},
      {198,       "T^{4}"},
      {199,       "T^{5}"},
      {200,   "T_{h}^{1}"},
      {201,   "T_{h}^{2}"},
      {202,   "T_{h}^{3}"},
      {203,   "T_{h}^{4}"},
      {204,   "T_{h}^{5}"},
      {205,   "T_{h}^{6}"},
      {206,   "T_{h}^{7}"},
      {207,       "O^{1}"},
      {208,       "O^{2}"},
      {209,       "O^{3}"},
      {210,       "O^{4}"},
      {211,       "O^{5}"},
      {212,       "O^{6}"},
      {213,       "O^{7}"},
      {214,       "O^{8}"},
      {215,   "T_{d}^{1}"},
      {216,   "T_{d}^{2}"},
      {217,   "T_{d}^{3}"},
      {218,   "T_{d}^{4}"},
      {219,   "T_{d}^{5}"},
      {220,   "T_{d}^{6}"},
      {221,   "O_{h}^{1}"},
      {222,   "O_{h}^{2}"},
      {223,   "O_{h}^{3}"},
      {224,   "O_{h}^{4}"},
      {225,   "O_{h}^{5}"},
      {226,   "O_{h}^{6}"},
      {227,   "O_{h}^{7}"},
      {228,   "O_{h}^{8}"},
      {229,   "O_{h}^{9}"},
      {230,  "O_{h}^{10}"}
  }; ///< spacegroup number -> Schoenflies notation

  // spacegroup number -> Hall notation (setting 1)
  static const std::map<int, std::string> sgn_hall_set1 = {
      {  1,              "P 1"},
      {  2,             "-P 1"},
      {  3,             "P 2y"},
      {  4,            "P 2yb"},
      {  5,             "C 2y"},
      {  6,            "P -2y"},
      {  7,           "P -2yc"},
      {  8,            "C -2y"},
      {  9,           "C -2yc"},
      { 10,            "-P 2y"},
      { 11,           "-P 2yb"},
      { 12,            "-C 2y"},
      { 13,           "-P 2yc"},
      { 14,          "-P 2ybc"},
      { 15,           "-C 2yc"},
      { 16,            "P 2 2"},
      { 17,           "P 2c 2"},
      { 18,          "P 2 2ab"},
      { 19,        "P 2ac 2ab"},
      { 20,           "C 2c 2"},
      { 21,            "C 2 2"},
      { 22,            "F 2 2"},
      { 23,            "I 2 2"},
      { 24,          "I 2b 2c"},
      { 25,           "P 2 -2"},
      { 26,          "P 2c -2"},
      { 27,          "P 2 -2c"},
      { 28,          "P 2 -2a"},
      { 29,        "P 2c -2ac"},
      { 30,         "P 2 -2bc"},
      { 31,         "P 2ac -2"},
      { 32,         "P 2 -2ab"},
      { 33,         "P 2c -2n"},
      { 34,          "P 2 -2n"},
      { 35,           "C 2 -2"},
      { 36,          "C 2c -2"},
      { 37,          "C 2 -2c"},
      { 38,           "A 2 -2"},
      { 39,          "A 2 -2c"},
      { 40,          "A 2 -2a"},
      { 41,         "A 2 -2ac"},
      { 42,           "F 2 -2"},
      { 43,          "F 2 -2d"},
      { 44,           "I 2 -2"},
      { 45,          "I 2 -2c"},
      { 46,          "I 2 -2a"},
      { 47,           "-P 2 2"},
      { 48,        "P 2 2 -1n"},
      { 49,          "-P 2 2c"},
      { 50,       "P 2 2 -1ab"},
      { 51,         "-P 2a 2a"},
      { 52,        "-P 2a 2bc"},
      { 53,         "-P 2ac 2"},
      { 54,        "-P 2a 2ac"},
      { 55,         "-P 2 2ab"},
      { 56,       "-P 2ab 2ac"},
      { 57,         "-P 2c 2b"},
      { 58,          "-P 2 2n"},
      { 59,     "P 2 2ab -1ab"},
      { 60,        "-P 2n 2ab"},
      { 61,       "-P 2ac 2ab"},
      { 62,        "-P 2ac 2n"},
      { 63,          "-C 2c 2"},
      { 64,         "-C 2bc 2"},
      { 65,           "-C 2 2"},
      { 66,          "-C 2 2c"},
      { 67,          "-C 2b 2"},
      { 68,       "C 2 2 -1bc"},
      { 69,           "-F 2 2"},
      { 70,        "F 2 2 -1d"},
      { 71,           "-I 2 2"},
      { 72,          "-I 2 2c"},
      { 73,         "-I 2b 2c"},
      { 74,          "-I 2b 2"},
      { 75,              "P 4"},
      { 76,             "P 4w"},
      { 77,             "P 4c"},
      { 78,            "P 4cw"},
      { 79,              "I 4"},
      { 80,            "I 4bw"},
      { 81,             "P -4"},
      { 82,             "I -4"},
      { 83,             "-P 4"},
      { 84,            "-P 4c"},
      { 85,       "P 4ab -1ab"},
      { 86,         "P 4n -1n"},
      { 87,             "-I 4"},
      { 88,       "I 4bw -1bw"},
      { 89,            "P 4 2"},
      { 90,        "P 4ab 2ab"},
      { 91,          "P 4w 2c"},
      { 92,       "P 4abw 2nw"},
      { 93,           "P 4c 2"},
      { 94,          "P 4n 2n"},
      { 95,         "P 4cw 2c"},
      { 96,       "P 4nw 2abw"},
      { 97,            "I 4 2"},
      { 98,        "I 4bw 2bw"},
      { 99,           "P 4 -2"},
      {100,         "P 4 -2ab"},
      {101,         "P 4c -2c"},
      {102,         "P 4n -2n"},
      {103,          "P 4 -2c"},
      {104,          "P 4 -2n"},
      {105,          "P 4c -2"},
      {106,        "P 4c -2ab"},
      {107,           "I 4 -2"},
      {108,          "I 4 -2c"},
      {109,         "I 4bw -2"},
      {110,        "I 4bw -2c"},
      {111,           "P -4 2"},
      {112,          "P -4 2c"},
      {113,         "P -4 2ab"},
      {114,          "P -4 2n"},
      {115,          "P -4 -2"},
      {116,         "P -4 -2c"},
      {117,        "P -4 -2ab"},
      {118,         "P -4 -2n"},
      {119,          "I -4 -2"},
      {120,         "I -4 -2c"},
      {121,           "I -4 2"},
      {122,         "I -4 2bw"},
      {123,           "-P 4 2"},
      {124,          "-P 4 2c"},
      {125,       "P 4 2 -1ab"},
      {126,        "P 4 2 -1n"},
      {127,         "-P 4 2ab"},
      {128,          "-P 4 2n"},
      {129,   "P 4ab 2ab -1ab"},
      {130,    "P 4ab 2n -1ab"},
      {131,          "-P 4c 2"},
      {132,         "-P 4c 2c"},
      {133,      "P 4n 2c -1n"},
      {134,       "P 4n 2 -1n"},
      {135,        "-P 4c 2ab"},
      {136,         "-P 4n 2n"},
      {137,      "P 4n 2n -1n"},
      {138,     "P 4n 2ab -1n"},
      {139,           "-I 4 2"},
      {140,          "-I 4 2c"},
      {141,   "I 4bw 2bw -1bw"},
      {142,   "I 4bw 2aw -1bw"},
      {143,              "P 3"},
      {144,             "P 31"},
      {145,             "P 32"},
      {146,             "P 3*"},
      {147,             "-P 3"},
      {148,            "-P 3*"},
      {149,            "P 3 2"},
      {150,          "P 3 2''"},
      {151,  "P 31 2c (0 0 1)"},
      {152,         "P 31 2''"},
      {153, "P 32 2c (0 0 -1)"},
      {154,         "P 32 2''"},
      {155,           "P 3* 2"},
      {156,         "P 3 -2''"},
      {157,           "P 3 -2"},
      {158,        "P 3 -2''c"},
      {159,          "P 3 -2c"},
      {160,          "P 3* -2"},
      {161,         "P 3* -2n"},
      {162,           "-P 3 2"},
      {163,          "-P 3 2c"},
      {164,         "-P 3 2''"},
      {165,        "-P 3 2''c"},
      {166,          "-P 3* 2"},
      {167,         "-P 3* 2n"},
      {168,              "P 6"},
      {169,             "P 61"},
      {170,             "P 65"},
      {171,             "P 62"},
      {172,             "P 64"},
      {173,             "P 6c"},
      {174,             "P -6"},
      {175,             "-P 6"},
      {176,            "-P 6c"},
      {177,            "P 6 2"},
      {178,  "P 61 2 (0 0 -1)"},
      {179,   "P 65 2 (0 0 1)"},
      {180,  "P 62 2c (0 0 1)"},
      {181, "P 64 2c (0 0 -1)"},
      {182,          "P 6c 2c"},
      {183,           "P 6 -2"},
      {184,          "P 6 -2c"},
      {185,          "P 6c -2"},
      {186,         "P 6c -2c"},
      {187,           "P -6 2"},
      {188,          "P -6c 2"},
      {189,          "P -6 -2"},
      {190,        "P -6c -2c"},
      {191,           "-P 6 2"},
      {192,          "-P 6 2c"},
      {193,          "-P 6c 2"},
      {194,         "-P 6c 2c"},
      {195,          "P 2 2 3"},
      {196,          "F 2 2 3"},
      {197,          "I 2 2 3"},
      {198,      "P 2ac 2ab 3"},
      {199,        "I 2b 2c 3"},
      {200,         "-P 2 2 3"},
      {201,      "P 2 2 3 -1n"},
      {202,         "-F 2 2 3"},
      {203,      "F 2 2 3 -1d"},
      {204,         "-I 2 2 3"},
      {205,     "-P 2ac 2ab 3"},
      {206,       "-I 2b 2c 3"},
      {207,          "P 4 2 3"},
      {208,         "P 4n 2 3"},
      {209,          "F 4 2 3"},
      {210,         "F 4d 2 3"},
      {211,          "I 4 2 3"},
      {212,     "P 4acd 2ab 3"},
      {213,      "P 4bd 2ab 3"},
      {214,       "I 4bd 2c 3"},
      {215,         "P -4 2 3"},
      {216,         "F -4 2 3"},
      {217,         "I -4 2 3"},
      {218,        "P -4n 2 3"},
      {219,        "F -4c 2 3"},
      {220,      "I -4bd 2c 3"},
      {221,         "-P 4 2 3"},
      {222,      "P 4 2 3 -1n"},
      {223,        "-P 4n 2 3"},
      {224,     "P 4n 2 3 -1n"},
      {225,         "-F 4 2 3"},
      {226,        "-F 4c 2 3"},
      {227,     "F 4d 2 3 -1d"},
      {228,    "F 4d 2 3 -1cd"},
      {229,         "-I 4 2 3"},
      {230,      "-I 4bd 2c 3"}
  }; ///< spacegroup number -> Hall notation (setting 1)

  // spacegroup number -> Hall notation (setting 2)
  static const std::map<int, std::string> sgn_hall_set2 = {
      { 48,    "-P 2ab 2bc"},
      { 50,     "-P 2ab 2b"},
      { 59,     "-P 2ab 2a"},
      { 68,     "-C 2b 2bc"},
      { 70,    "-F 2uv 2vw"},
      { 85,         "-P 4a"},
      { 86,        "-P 4bc"},
      { 88,        "-I 4ad"},
      {125,      "-P 4a 2b"},
      {126,     "-P 4a 2bc"},
      {129,      "-P 4a 2a"},
      {130,     "-P 4a 2ac"},
      {133,     "-P 4ac 2b"},
      {134,    "-P 4ac 2bc"},
      {137,     "-P 4ac 2a"},
      {138,    "-P 4ac 2ac"},
      {141,      "-I 4bd 2"},
      {142,     "-I 4bd 2c"},
      {146,           "R 3"},
      {148,          "-R 3"},
      {155,       "R 3 2''"},
      {160,      "R 3 -2''"},
      {161,     "R 3 -2''c"},
      {166,      "-R 3 2''"},
      {167,     "-R 3 2''c"},
      {201,  "-P 2ab 2bc 3"},
      {203,  "-F 2uv 2vw 3"},
      {222,   "-P 4a 2bc 3"},
      {224,  "-P 4bc 2bc 3"},
      {227,  "-F 4vw 2vw 3"},
      {228, "-F 4cvw 2vw 3"}
  }; ///< spacegroup number -> Hall notation (setting 2)

} // namespace symmetry_lookup

// ***************************************************************************
// GetSpaceGroupName
// ***************************************************************************
string GetSpaceGroupName(const int spacegroupnumber, const string& directory) {
  if (spacegroupnumber < 1 || spacegroupnumber > 230) { // DX20190708 - for xerror
    stringstream message; // DX20190708 - for xerror
    message << "routine: space group specified invalid (1-230): "; // DX20190708 - for xerror
    message << spacegroupnumber << " [dir=" << directory << "]." << endl; // DX20190708 - for xerror
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_); // DX20190708 - for xerror
  }

  // HE20240530
  if (symmetry_lookup::sgn_sg.count(spacegroupnumber)) {
    return symmetry_lookup::sgn_sg.at(spacegroupnumber); // .at instead of [] needed as map is constant
  }
  return "";
}

// ***************************************************************************
// GetSpaceGroupNumber
// ***************************************************************************
int GetSpaceGroupNumber(const string& spacegroupsymbol, const string& directory) {
  // DX20190708
  stringstream message;
  if (spacegroupsymbol[0] != 'P' && spacegroupsymbol[0] != 'I' && spacegroupsymbol[0] != 'F' && spacegroupsymbol[0] != 'R' && spacegroupsymbol[0] != 'C' && spacegroupsymbol[0] != 'A') {
    message << "routine: space group specified invalid (lattice centering not identified: P,I,F,R,C,A): ";
    message << "input symbol=" << spacegroupsymbol << " [dir=" << directory << "]." << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }

  // HE20240530
  if (symmetry_lookup::sg_sgn.count(spacegroupsymbol)) {
    return symmetry_lookup::sg_sgn.at(spacegroupsymbol); // .at instead of [] needed as map is constant
  } else {
    message << "routine: space group specified invalid; perhaps non-ITC setting: ";
    message << "space group symbol=" << spacegroupsymbol << " [dir=" << directory << "].";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
}

// ***************************************************************************
// GetSpaceGroupSchoenflies
// ***************************************************************************
string GetSpaceGroupSchoenflies(const int spacegroupnumber, const string& directory) {
  if (spacegroupnumber < 1 || spacegroupnumber > 230) { // DX20190708 - for xerror
    stringstream message; // DX20190708 - for xerror
    message << "routine: space group specified invalid (1-230): "; // DX20190708 - for xerror
    message << spacegroupnumber << " [dir=" << directory << "]." << endl; // DX20190708 - for xerror
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_); // DX20190708 - for xerror
  }
  // HE20240530
  if (symmetry_lookup::sgn_sch.count(spacegroupnumber)) {
    return symmetry_lookup::sgn_sch.at(spacegroupnumber); // .at instead of [] needed as map is constant
  }
  // No space group found
  return "";
}

// ***************************************************************************
// GetSpaceGroupHall
// ***************************************************************************
string GetSpaceGroupHall(const int spacegroupnumber, int setting, const string& directory) {
  // DX - Hall distinguishes space group setting.  This table assumes the first
  //       setting that appears in the ITC.
  //       For more settings, they need to be hard-coded here.
  stringstream message; // DX20190708 - for xerror
  if (spacegroupnumber < 1 || spacegroupnumber > 230) { // DX20190708 - for xerror
    message << "routine: space group specified invalid (1-230): "; // DX20190708 - for xerror
    message << spacegroupnumber << " [dir=" << directory << "]." << endl; // DX20190708 - for xerror
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_); // DX20190708 - for xerror
  }
  // OK
  if (setting == SG_SETTING_ANRL) {
    setting = anrl::getANRLSettingChoice(spacegroupnumber);
  } // DX20210420
  else if (setting == 0) { // signals default //DX20180807
    // if RHL, AFLOW prefers hexagonal setting (i.e., setting=2)
    if (spacegroupnumber == 146 || spacegroupnumber == 148 || spacegroupnumber == 155 || spacegroupnumber == 160 || spacegroupnumber == 161 || spacegroupnumber == 166 || spacegroupnumber == 167) {
      setting = 2;
    }
    // else setting==SG_SETTING_1
    else {
      setting = SG_SETTING_1;
    }
  }
  if (setting < 1 || setting > 2) {
    message << "routine: setting choice is invalid (1 or 2 only): " << setting << " [dir=" << directory << "].";
    // DX20190708 - for xerror
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_); // DX20190708 - for xerror
  }

  // HE20240530
  if (setting == SG_SETTING_2) {
    if (symmetry_lookup::sgn_hall_set2.count(spacegroupnumber)) {
      return symmetry_lookup::sgn_hall_set2.at(spacegroupnumber);
    }
  }
  if (symmetry_lookup::sgn_hall_set1.count(spacegroupnumber)) {
    return symmetry_lookup::sgn_hall_set1.at(spacegroupnumber);
  }
  // No space group found
  return "";
}

// ***************************************************************************
// GetLaueLabel
// ***************************************************************************
string GetLaueLabel(string& point_group) {
  string laue;
  // DX+ME20190708 - changed subsequent "if" to "else if" -> efficiency
  //  -1
  if (point_group == "1" || point_group == "-1") {
    laue = "-1";
  }
  // 2/m
  else if (point_group == "2" || point_group == "m" || point_group == "2/m") {
    laue = "2/m";
  }
  // mmm
  else if (point_group == "222" || point_group == "mm2" || point_group == "mmm") {
    laue = "mmm";
  }
  // 4/m
  else if (point_group == "4" || point_group == "-4" || point_group == "4/m") {
    laue = "4/m";
  }
  // 4/mmmm
  else if (point_group == "422" || point_group == "4mm" || point_group == "-42m" || point_group == "-4m2" || point_group == "4/mmm") {
    laue = "4/mmm";
  }
  // -3
  else if (point_group == "3" || point_group == "-3") {
    laue = "-3";
  }
  // -3m
  else if (point_group == "312" || point_group == "321" || point_group == "32" || point_group == "31m" || point_group == "3m1" || point_group == "-31m" || point_group == "-3m1" || point_group == "-3m" ||
           point_group == "3m") {
    laue = "-3m";
  }
  // 6/m
  else if (point_group == "6" || point_group == "-6" || point_group == "6/m") {
    laue = "6/m";
  }
  // 6/mmm
  else if (point_group == "622" || point_group == "6mm" || point_group == "-6m2" || point_group == "-62m" || point_group == "6/mmm") {
    laue = "6/mmm";
  }
  // m-3
  else if (point_group == "23" || point_group == "m-3") {
    laue = "m-3";
  }
  // m-3m
  else if (point_group == "432" || point_group == "-43m" || point_group == "m-3m") {
    laue = "m-3m";
  }
  return laue;
}

// ***************************************************************************
// GetSpaceGroupLabel
// ***************************************************************************
string GetSpaceGroupLabel(int spacegroupnumber) {
  string spacegrouplabel;
  spacegrouplabel = "#" + aurostd::utype2string(spacegroupnumber);
  return spacegrouplabel;
}

// **************************************************************************
// Function MetricTensor
// **************************************************************************
// this function returns the metric tensor
// Corey Oses
xmatrix<double> MetricTensor(const xstructure& a) {
  return MetricTensor(a.lattice, a.scale);
}

xmatrix<double> MetricTensor(const xmatrix<double>& lattice, double scale) {
  if (lattice.rows != lattice.cols) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Dimension mismatch, should be square lattice matrix.", _VALUE_ILLEGAL_);
  }
  const xmatrix<double> metric_tensor(lattice.rows, lattice.cols);
  for (int i = lattice.lrows; i <= lattice.urows; i++) { // CO20190520
    for (int j = lattice.lcols; j <= lattice.ucols; j++) { // CO20190520
      metric_tensor(i, j) = aurostd::scalar_product(lattice(i), lattice(j));
    }
  }
  return scale * metric_tensor;
}

// **************************************************************************
// Function ReciprocalLattice
// **************************************************************************
// this function returns the reciprocal lattice
// from the original one... it takes the scale !
// Stefano Curtarolo
xmatrix<double> ReciprocalLattice(const xstructure& a) {
  return ReciprocalLattice(a.lattice, a.scale);
}

xmatrix<double> ReciprocalLattice(const xmatrix<double>& rlattice, double scale) { // AFLOW_FUNCTION_IMPLEMENTATION
  const xvector<double> a1(3);
  const xvector<double> a2(3);
  const xvector<double> a3(3);
  xvector<double> b1(3);
  xvector<double> b2(3);
  xvector<double> b3(3);
  const xmatrix<double> klattice(3, 3); // kvectors are RAWS
  double norm;
  a1[1] = rlattice[1][1];
  a1[2] = rlattice[1][2];
  a1[3] = rlattice[1][3];
  a2[1] = rlattice[2][1];
  a2[2] = rlattice[2][2];
  a2[3] = rlattice[2][3];
  a3[1] = rlattice[3][1];
  a3[2] = rlattice[3][2];
  a3[3] = rlattice[3][3];
  norm = 2.0 * pi / (det(rlattice) * scale);
  b1 = norm * vector_product(a2, a3);
  b2 = norm * vector_product(a3, a1);
  b3 = norm * vector_product(a1, a2);
  klattice[1][1] = b1[1];
  klattice[1][2] = b1[2];
  klattice[1][3] = b1[3];
  klattice[2][1] = b2[1];
  klattice[2][2] = b2[2];
  klattice[2][3] = b2[3];
  klattice[3][1] = b3[1];
  klattice[3][2] = b3[2];
  klattice[3][3] = b3[3];
  return klattice;
}

// xmatrix<double> ReciprocalLattice(const xmatrix<double>& rlattice) {        // AFLOW_FUNCTION_IMPLEMENTATION
//   return ReciprocalLattice(rlattice,1.0);
// }

// **************************************************************************
// Function KPPRA
// **************************************************************************
// This function calculates k1,k2,k3 starting from real lattice and NK total
// the function does not normalize with number of atoms so the calculation
// must be done somewhere else
string KPPRA(int& k1, int& k2, int& k3, const xmatrix<double>& rlattice, const int& NK) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xmatrix<double> klattice(3, 3);
  klattice = ReciprocalLattice(rlattice, 1.0);
  xvector<double> b1(3);
  xvector<double> b2(3);
  xvector<double> b3(3);
  xvector<double> db1(3);
  xvector<double> db2(3);
  xvector<double> db3(3);
  double nb1;
  double nb2;
  double nb3;
  k1 = 1;
  k2 = 1;
  k3 = 1;
  int kk1;
  int kk2;
  int kk3;
  int kk;
  b1 = klattice(1);
  nb1 = aurostd::modulus(b1);
  b2 = klattice(2);
  nb2 = aurostd::modulus(b2);
  b3 = klattice(3);
  nb3 = aurostd::modulus(b3);
  if (LDEBUG) {
    aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  }
  if (LDEBUG) {
    aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  }
  if (NK > 1) {
    bool found = false;
    double dkdelta;
    double dk;
    dkdelta = 0.999;
    dk = aurostd::min(nb1, nb2, nb3);
    kk1 = 0;
    kk2 = 0;
    kk3 = 0;
    kk = 0;
    int iverbose = 0;
    while (!found) {
      kk++;
      if (dk <= 1e-5) {
        k1 = 1;
        k2 = 1;
        k3 = 1;
        aus << "00000  MESSAGE KPOINTS KPPRA minimum not found k=[" << k1 << "," << k2 << "," << k3 << "]=" << k1 * k2 * k3 << endl;
      }
      kk1 = (int) floor(nb1 / dk);
      db1 = b1 / ((double) kk1);
      kk2 = (int) floor(nb2 / dk);
      db2 = b2 / ((double) kk2);
      kk3 = (int) floor(nb3 / dk);
      db3 = b3 / ((double) kk3);
      if (kk1 + kk2 + kk3 > iverbose) {
        //  if(!mod(kk,50) || kk1*kk2*kk3>=NK)
        aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "]=" << kk1 * kk2 * kk3 << " =[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << ","
            << aurostd::modulus(db3) << "]   dk=" << dk << endl;
        iverbose = kk1 + kk2 + kk3;
      }
      if (kk1 * kk2 * kk3 >= NK) {
        k1 = kk1;
        k2 = kk2;
        k3 = kk3;
        found = true;
      }
      dk = dk * dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1 = 1;
    k2 = 1;
    k3 = 1;
  }
  db1 = b1 / ((double) k1);
  db2 = b2 / ((double) k2);
  db3 = b3 / ((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1 * k2 * k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   " << endl;
  return aus.str();
}

string KPPRA_LAT(int& k1, int& k2, int& k3, const xmatrix<double>& rlattice, const int& NK) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xvector<double> kdata(6);
  xvector<double> kdatagrid(6);
  xmatrix<double> klattice(3, 3);
  klattice = ReciprocalLattice(rlattice, 1.0);
  klattice = GetStandardPrimitive(klattice);
  kdata = Getabc_angles(klattice, DEGREES);

  xvector<double> b1(3);
  xvector<double> b2(3);
  xvector<double> b3(3);
  xvector<double> db1(3);
  xvector<double> db2(3);
  xvector<double> db3(3);
  double nb1;
  double nb2;
  double nb3;
  k1 = 1;
  k2 = 1;
  k3 = 1;
  int kk1;
  int kk2;
  int kk3;
  int kk;
  b1 = klattice(1);
  nb1 = aurostd::modulus(b1);
  b2 = klattice(2);
  nb2 = aurostd::modulus(b2);
  b3 = klattice(3);
  nb3 = aurostd::modulus(b3);
  if (LDEBUG) {
    aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  }
  if (LDEBUG) {
    aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  }
  if (NK > 1) {
    bool found = false;
    double dkdelta;
    double dk;
    dkdelta = 0.999;
    dk = aurostd::min(nb1, nb2, nb3);
    kk1 = 0;
    kk2 = 0;
    kk3 = 0;
    kk = 0;
    int iverbose = 0;
    while (!found) {
      kk++;
      if (dk <= 1e-5) {
        k1 = 1;
        k2 = 1;
        k3 = 1;
        aus << "00000  MESSAGE KPOINTS KPPRA minimum not found k=[" << k1 << "," << k2 << "," << k3 << "]=" << k1 * k2 * k3 << endl;
      }
      kk1 = (int) floor(nb1 / dk);
      db1 = b1 / ((double) kk1);
      kk2 = (int) floor(nb2 / dk);
      db2 = b2 / ((double) kk2);
      kk3 = (int) floor(nb3 / dk);
      db3 = b3 / ((double) kk3);
      kdatagrid = kdata;
      kdatagrid[1] = 1.0;
      kdatagrid[2] = 1.0 * aurostd::modulus(db2) / aurostd::modulus(db1);
      kdatagrid[3] = 1.0 * aurostd::modulus(db3) / aurostd::modulus(db1);

      if (kk1 + kk2 + kk3 > iverbose) {
        //  if(!mod(kk,50) || kk1*kk2*kk3>=NK)
        aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "," << GetLatticeType(kdatagrid) << "]=" << kk1 * kk2 * kk3 << " = [" << aurostd::modulus(db1) << ","
            << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   dk=" << dk << " " << endl;
        iverbose = kk1 + kk2 + kk3;
      }
      if (kk1 * kk2 * kk3 >= NK) {
        k1 = kk1;
        k2 = kk2;
        k3 = kk3;
        found = true;
      }
      dk = dk * dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1 = 1;
    k2 = 1;
    k3 = 1;
  }
  db1 = b1 / ((double) k1);
  db2 = b2 / ((double) k2);
  db3 = b3 / ((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1 * k2 * k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   " << endl;
  return aus.str();
}

string KPPRA(xstructure& str, const int& _NK) {
  //  cerr << "KPPRA" << endl;
  int NK = 1;
  NK = (int) ((double) _NK / str.atoms.size() + 0.5);
  if (NK < 1) {
    NK = 1;
  }
  int k1 = 1;
  int k2 = 1;
  int k3 = 1;
  xmatrix<double> rlattice = str.lattice;
  rlattice = str.scale * rlattice;
  // string stringKPPRA=KPPRA_LAT(k1,k2,k3,rlattice,NK);
  string stringKPPRA = KPPRA(k1, k2, k3, rlattice, NK);
  str.kpoints_k1 = k1;
  str.kpoints_k2 = k2;
  str.kpoints_k3 = k3;
  str.kpoints_kmax = max(str.kpoints_k1, str.kpoints_k2, str.kpoints_k3);
  str.kpoints_kppra = str.kpoints_k1 * str.kpoints_k2 * str.kpoints_k3 * str.atoms.size();
  return stringKPPRA;
}

// **************************************************************************
// Function KPPRA_DELTA
// **************************************************************************
// This function calculates k1,k2,k3 starting from real lattice and DK
string KPPRA_DELTA(int& k1, int& k2, int& k3, const xmatrix<double>& rlattice, const double& DK) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream aus("");
  aus.precision(5);
  xmatrix<double> klattice(3, 3);
  klattice = ReciprocalLattice(rlattice, 1.0);
  xvector<double> b1(3);
  xvector<double> b2(3);
  xvector<double> b3(3);
  xvector<double> db1(3);
  xvector<double> db2(3);
  xvector<double> db3(3);
  double nb1;
  double nb2;
  double nb3;
  k1 = 1;
  k2 = 1;
  k3 = 1;
  int kk1;
  int kk2;
  int kk3;
  int kk;
  b1 = klattice(1);
  nb1 = aurostd::modulus(b1);
  b2 = klattice(2);
  nb2 = aurostd::modulus(b2);
  b3 = klattice(3);
  nb3 = aurostd::modulus(b3);
  if (LDEBUG) {
    aus << "KPPRA LDEBUG:  " << endl << rlattice << endl << endl << klattice << endl << b1 << endl << b2 << endl << b3 << endl << nb1 << endl << nb2 << endl << nb3 << endl;
  }
  if (LDEBUG) {
    aus << "KPPRA LDEBUG:  " << nb1 << " " << nb2 << " " << nb3 << " " << endl;
  }
  if (DK > 1.0e-6) {
    bool found = false;
    double dkdelta;
    double dk;
    dkdelta = 0.9999;
    dk = aurostd::min(nb1, nb2, nb3);
    kk1 = 0;
    kk2 = 0;
    kk3 = 0;
    kk = 0;
    int iverbose = 0;
    while (!found) {
      kk++;
      kk1 = (int) floor(nb1 / dk);
      db1 = b1 / ((double) kk1);
      kk2 = (int) floor(nb2 / dk);
      db2 = b2 / ((double) kk2);
      kk3 = (int) floor(nb3 / dk);
      db3 = b3 / ((double) kk3);
      if (kk1 + kk2 + kk3 > iverbose) {
        //  if(!mod(kk,50) || kk1*kk2*kk3>=DK)
        aus << "00000  MESSAGE KPOINTS KPPRA minimizing k=[" << kk1 << "," << kk2 << "," << kk3 << "]=" << kk1 * kk2 * kk3 << " =[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << ","
            << aurostd::modulus(db3) << "]   dk=" << dk << endl;
        iverbose = kk1 + kk2 + kk3;
      }
      if ((aurostd::modulus(db1) < DK) && (aurostd::modulus(db2) < DK) && (aurostd::modulus(db3) < DK)) {
        found = true;
        k1 = kk1;
        k2 = kk2;
        k3 = kk3;
      }
      dk = dk * dkdelta;
    }
  } else { // force 1 1 1 for Gamma
    k1 = 1;
    k2 = 1;
    k3 = 1;
  }
  db1 = b1 / ((double) k1);
  db2 = b2 / ((double) k2);
  db3 = b3 / ((double) k3);
  aus << "00000  MESSAGE KPOINTS KPPRA routine [" << k1 << "," << k2 << "," << k3 << "]=" << k1 * k2 * k3 << "=[" << aurostd::modulus(db1) << "," << aurostd::modulus(db2) << "," << aurostd::modulus(db3) << "]   " << endl;
  return aus.str();
}

string KPPRA_DELTA(xstructure& str, const double& DK) {
  //  cerr << "KPPRA_DELTA" << endl;
  int k1 = 1;
  int k2 = 1;
  int k3 = 1;
  xmatrix<double> rlattice = str.lattice;
  rlattice = str.scale * rlattice;
  string stringKPPRA = KPPRA_DELTA(k1, k2, k3, rlattice, DK);
  str.kpoints_k1 = k1;
  str.kpoints_k2 = k2;
  str.kpoints_k3 = k3;
  str.kpoints_kmax = max(str.kpoints_k1, str.kpoints_k2, str.kpoints_k3);
  str.kpoints_kppra = str.kpoints_k1 * str.kpoints_k2 * str.kpoints_k3 * str.atoms.size();
  return stringKPPRA;
}

// **************************************************************************
// Function GetNBAND
// **************************************************************************
// returns estimated version of NBANDS starting from
// electrons, ions, spin and ispin
int GetNBANDS(int electrons, int nions, int spineach, bool ispin, int NPAR) {
  double out = 0.0;
  out = max(ceil((electrons + 4.0) / 1.75) + max(nions / 1.75, 6.0), ceil(0.80 * electrons)); // from VASP
  if (ispin) {
    out += (nions * spineach + 1) / 2;
  }
  //  out*=1.2;  // safety from vasp
  out *= 1.3; // safety more
  out = out * 1.1 + 5; // Thu Jun 11 12:08:42 EDT 2009 // METAL PROJECT
  out *= 1.075; // Tue Oct 13 07:59:43 EDT 2009 // ICSD PROJECT
  out += 5; // Sun Nov  1 10:41:20 EDT 2009 // ICSD PROJECT ORC
  out *= 1.03; // Tue Feb 26 15:15:36 EST 2013 // HELPS dielectric CALS
  out *= 1.05; // Mon Apr 23 13:40:02 EST 2018 // HELPS SCAN
  if (nions < 100) {
    out *= std::pow((double) nions, 0.025);
    // rescale so for big numbers of ions you get extra bands // Wed Jun 23 12:29:01 EDT 2010
  } else {
    out *= std::pow((double) nions, 0.06); // ME20191028 - prior scaling factor not sufficient for supercells
  }
  int nbands = (int) ceil(out);
  // CO20210315 START - adjust for NPAR
  if (NPAR > 0) {
    const int increment = 1; //(increase?+1:-1);
    while ((nbands % NPAR) != 0) {
      nbands += increment;
    }
  }
  // CO20210315 END - adjust for NPAR
  return nbands;
}

// **************************************************************************
// Function GetZVAL from *CAR
// ***************************************************************************
double GetZVAL(const stringstream& sss, vector<double>& vZVAL) {
  xPOTCAR potcar;
  potcar.GetProperties(sss);
  vZVAL.clear();
  for (size_t i = 0; i < potcar.vZVAL.size(); i++) {
    vZVAL.push_back(potcar.vZVAL[i]);
  }
  return potcar.ZVAL_sum;
}

double GetZVAL(const _xvasp& xvasp, vector<double>& vZVAL) {
  return GetZVAL(xvasp.POTCAR, vZVAL);
}

double GetZVAL(const string& directory, vector<double>& vZVAL) {
  stringstream sss("");
  const vector<string> vfile{"POTCAR", "OUTCAR", "POTCAR.relax1", "POTCAR.relax2", "POTCAR.static", "POTCAR.bands", "OUTCAR.relax1", "OUTCAR.relax2", "OUTCAR.static", "OUTCAR.bands"};
  for (size_t i = 0; i < vfile.size() && sss.str().empty(); i++) {
    if (sss.str().empty() && aurostd::FileExist(directory + "/" + vfile[i])) {
      aurostd::file2stringstream(directory + "/" + vfile[i], sss);
    }
    if (sss.str().empty() && aurostd::CompressFileExist(directory + "/" + vfile[i])) {
      aurostd::compressfile2stringstream(directory + "/" + vfile[i], sss);
    }
  } //  cerr << sss.str() << endl;
  return GetZVAL(sss, vZVAL);
}

double GetCellAtomZVAL(const stringstream& sss, vector<double>& vZVAL, const stringstream& sstr, vector<double>& sZVAL, string mode) {
  vZVAL.clear();
  sZVAL.clear();
  GetZVAL(sss, vZVAL);
  stringstream aus(sstr.str());
  aus << sstr.str();
  xstructure xstr(aus, IOAFLOW_AUTO);
  if (mode == "CELL" || mode.empty()) {
    const double CellZVAL = xstr.GetZVAL(vZVAL);
    for (size_t i = 0; i < vZVAL.size(); i++) {
      sZVAL.push_back(double(vZVAL[i] * double(xstr.num_each_type.at(i))));
    }
    return CellZVAL;
  }
  if (mode == "ATOM") {
    const double CellZVAL = xstr.GetZVAL(vZVAL) / double(xstr.atoms.size());
    for (size_t i = 0; i < vZVAL.size(); i++) {
      sZVAL.push_back(double(vZVAL[i] * double(xstr.num_each_type.at(i))) / double(xstr.atoms.size()));
    }
    return CellZVAL;
  }
  return 0.0;
}

double GetCellAtomZVAL(const string& directory, vector<double>& vZVAL, vector<double>& sZVAL, string mode) {
  // from directory POT/POS returns total ZVAL cell, vZVAL and sZVAL
  // search for data
  stringstream sss("");
  const vector<string> vfile1{"POTCAR", "OUTCAR", "POTCAR.relax1", "POTCAR.relax2", "POTCAR.static", "POTCAR.bands", "OUTCAR.relax1", "OUTCAR.relax2", "OUTCAR.static", "OUTCAR.bands"};
  for (size_t i = 0; i < vfile1.size() && sss.str().empty(); i++) {
    if (sss.str().empty() && aurostd::FileExist(directory + "/" + vfile1[i])) {
      aurostd::file2stringstream(directory + "/" + vfile1[i], sss);
    }
    if (sss.str().empty() && aurostd::CompressFileExist(directory + "/" + vfile1[i])) {
      aurostd::compressfile2stringstream(directory + "/" + vfile1[i], sss);
    }
  } //   cerr << sss.str() << endl;

  // search for xstructure
  stringstream sstr("");
  const vector<string> vfile2{"POSCAR", "CONTCAR", "POSCAR.relax1", "POSCAR.relax2", "POSCAR.static", "POSCAR.bands", "CONTCAR.relax1", "CONTCAR.relax2", "CONTCAR.static", "CONTCAR.bands"};
  for (size_t i = 0; i < vfile2.size() && sstr.str().empty(); i++) {
    //    cerr << "file=" << directory+"/"+vfile.at(i) << endl;
    if (sstr.str().empty() && aurostd::FileExist(directory + "/" + vfile2[i])) {
      aurostd::file2stringstream(directory + "/" + vfile2[i], sstr);
    }
    if (sstr.str().empty() && aurostd::CompressFileExist(directory + "/" + vfile2[i])) {
      aurostd::compressfile2stringstream(directory + "/" + vfile2[i], sstr);
    }
  } //  cerr << sstr.str() << endl;
  // done
  return GetCellAtomZVAL(sss, vZVAL, sstr, sZVAL, mode);
}

// ***************************************************************************
// Function GetZVAL
// ***************************************************************************
// Given the ZVAL of each species, it returns total ZVAL of cell
double xstructure::GetZVAL(const vector<double>& vZVAL) {
  stringstream message;
  if (num_each_type.size() != vZVAL.size()) {
    message << "num_each_type.size()=" << num_each_type.size() << endl;
    message << "vZVAL.size()=" << vZVAL.size() << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
  }
  double CellZVAL = 0.0;
  for (size_t i = 0; i < vZVAL.size(); i++) {
    CellZVAL += double(vZVAL[i] * num_each_type.at(i));
  }
  return CellZVAL;
}

// **************************************************************************
// Function GetPOMASS from *CAR
// ***************************************************************************
double GetPOMASS(const stringstream& sss, vector<double>& vPOMASS) {
  xPOTCAR potcar;
  potcar.GetProperties(sss);
  vPOMASS.clear();
  for (size_t i = 0; i < potcar.vPOMASS.size(); i++) {
    vPOMASS.push_back(potcar.vPOMASS[i]);
  }
  return potcar.POMASS_sum;
}

double GetPOMASS(const _xvasp& xvasp, vector<double>& vPOMASS) {
  return GetPOMASS(xvasp.POTCAR, vPOMASS);
}

double GetPOMASS(const string& directory, vector<double>& vPOMASS) {
  stringstream sss("");
  const vector<string> vfile{"POTCAR", "OUTCAR", "POTCAR.relax1", "POTCAR.relax2", "POTCAR.static", "POTCAR.bands", "OUTCAR.relax1", "OUTCAR.relax2", "OUTCAR.static", "OUTCAR.bands"};
  for (size_t i = 0; i < vfile.size() && sss.str().empty(); i++) {
    if (sss.str().empty() && aurostd::FileExist(directory + "/" + vfile[i])) {
      aurostd::file2stringstream(directory + "/" + vfile[i], sss);
    }
    if (sss.str().empty() && aurostd::CompressFileExist(directory + "/" + vfile[i])) {
      aurostd::compressfile2stringstream(directory + "/" + vfile[i], sss);
    }
  } //  cerr << sss.str() << endl;
  return GetPOMASS(sss, vPOMASS);
}

double GetCellAtomPOMASS(const stringstream& sss, vector<double>& vPOMASS, const stringstream& sstr, vector<double>& sPOMASS, string mode) {
  vPOMASS.clear();
  sPOMASS.clear();
  GetPOMASS(sss, vPOMASS);
  stringstream aus(sstr.str());
  aus << sstr.str();
  xstructure xstr(aus, IOAFLOW_AUTO);
  if (mode == "CELL" || mode.empty()) {
    const double CellPOMASS = xstr.GetPOMASS(vPOMASS);
    for (size_t i = 0; i < vPOMASS.size(); i++) {
      sPOMASS.push_back(double(vPOMASS[i] * double(xstr.num_each_type.at(i))));
    }
    return CellPOMASS;
  }
  if (mode == "ATOM") {
    const double CellPOMASS = xstr.GetPOMASS(vPOMASS) / double(xstr.atoms.size());
    for (size_t i = 0; i < vPOMASS.size(); i++) {
      sPOMASS.push_back(double(vPOMASS[i] * double(xstr.num_each_type.at(i))) / double(xstr.atoms.size()));
    }
    return CellPOMASS;
  }
  return 0.0;
}

double GetCellAtomPOMASS(const string& directory, vector<double>& vPOMASS, vector<double>& sPOMASS, string mode) {
  // from directory POT/POS returns total POMASS cell, vPOMASS and sPOMASS
  const vector<string> vfile1{
      "POTCAR", "OUTCAR", "POTCAR.relax1", "POTCAR.relax2", "POTCAR.static", "POTCAR.bands", "OUTCAR.relax1", "OUTCAR.relax2", "OUTCAR.static", "OUTCAR.bands",
  };
  // search for data
  stringstream sss("");

  for (size_t i = 0; i < vfile1.size() && sss.str().empty(); i++) {
    if (sss.str().empty() && aurostd::FileExist(directory + "/" + vfile1[i])) {
      aurostd::file2stringstream(directory + "/" + vfile1[i], sss);
    }
    if (sss.str().empty() && aurostd::CompressFileExist(directory + "/" + vfile1[i])) {
      aurostd::compressfile2stringstream(directory + "/" + vfile1[i], sss);
    }
  } //  cerr << sss.str() << endl;

  // search for xstructure
  stringstream sstr("");
  const vector<string> vfile2{"POSCAR", "CONTCAR", "POSCAR.relax1", "POSCAR.relax2", "POSCAR.static", "POSCAR.bands", "CONTCAR.relax1", "CONTCAR.relax2", "CONTCAR.static", "CONTCAR.bands"};

  for (size_t i = 0; i < vfile2.size() && sstr.str().empty(); i++) {
    if (sstr.str().empty() && aurostd::FileExist(directory + "/" + vfile2[i])) {
      aurostd::file2stringstream(directory + "/" + vfile2[i], sstr);
    }
    if (sstr.str().empty() && aurostd::CompressFileExist(directory + "/" + vfile2[i])) {
      aurostd::compressfile2stringstream(directory + "/" + vfile2[i], sstr);
    }
  } //  cerr << sstr.str() << endl;

  // done
  return GetCellAtomPOMASS(sss, vPOMASS, sstr, sPOMASS, mode);
}

// ***************************************************************************
// Function GetPOMASS
// ***************************************************************************
// Given the POMASS of each species, it returns total POMASS of cell
double xstructure::GetPOMASS(const vector<double>& vPOMASS) {
  stringstream message;
  if (num_each_type.size() != vPOMASS.size()) {
    message << "num_each_type.size()=" << num_each_type.size() << endl;
    message << "vPOMASS.size()=" << vPOMASS.size() << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
  }
  double CellPOMASS = 0.0;
  for (size_t i = 0; i < vPOMASS.size(); i++) {
    CellPOMASS += double(vPOMASS[i] * num_each_type.at(i));
  }
  return CellPOMASS;
}

// ***************************************************************************
// AtomEnvironment Class - DX20191122
// ***************************************************************************
// ---------------------------------------------------------------------------
// AtomEnvironment (constructor)
AtomEnvironment::AtomEnvironment() {
  free();
}

// ---------------------------------------------------------------------------
// AtomEnvironment::free
void AtomEnvironment::free() {
  mode = 0;
  element_center = "";
  type_center = 0;
  num_types = 0;
  num_neighbors = 0;
  elements_neighbor.clear();
  types_neighbor.clear();
  distances_neighbor.clear();
  coordinations_neighbor.clear();
  coordinates_neighbor.clear();
  facets.clear();
  facet_area.clear();
  area = 0;
  volume = 0;
  has_hull = false;
  facet_order.clear();
  facet_order.resize(8, 0);
}

// ---------------------------------------------------------------------------
// AtomEnvironment (destructor)
AtomEnvironment::~AtomEnvironment() {
  free();
}

// ---------------------------------------------------------------------------
// AtomEnvironment (copy constructor)
AtomEnvironment::AtomEnvironment(const AtomEnvironment& b) {
  copy(b);
}

// ---------------------------------------------------------------------------
// AtomEnvironment::copy
void AtomEnvironment::copy(const AtomEnvironment& b) {
  mode = b.mode;
  element_center = b.element_center;
  type_center = b.type_center;
  num_neighbors = b.num_neighbors;
  num_types = b.num_types;
  elements_neighbor = b.elements_neighbor;
  types_neighbor = b.types_neighbor;
  distances_neighbor = b.distances_neighbor;
  coordinations_neighbor = b.coordinations_neighbor;
  coordinates_neighbor = b.coordinates_neighbor;
  facets = b.facets;
  facet_order = b.facet_order;
  facet_area = b.facet_area;
  area = b.area;
  volume = b.volume;
  has_hull = b.has_hull;
}

// ---------------------------------------------------------------------------
// AtomEnvironment::operator=
const AtomEnvironment& AtomEnvironment::operator=(const AtomEnvironment& b) {
  if (this != &b) {
    copy(b);
  }
  return *this;
}

// ---------------------------------------------------------------------------
// AtomEnvironment::operator<<
ostream& operator<<(ostream& cout, const AtomEnvironment& AtomEnvironment) {
  cout << AtomEnvironment.toJSON(false).toString();

  return cout;
}

// ***************************************************************************
// AtomEnvironment::constructAtomEnvironmentHull() - HE20210408
// ***************************************************************************

/// @brief constructed a convex hull around the atomic environment
void AtomEnvironment::constructAtomEnvironmentHull() {
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (has_hull) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AE hull is already set" << endl;
    }
    return;
  }

  vector<xvector<double>> points;
  for (uint t = 0; t < num_neighbors; t++) {
    points.push_back(index2Point(t));
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " create AE hull around " << num_neighbors << " atoms" << endl;
  }
  xoption hull_options;
  hull_options.flag("CHULL::FULL_HULL", true);
  hull_options.flag("CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS", true);
  hull_options.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS", true);
  hull_options.flag("CHULL::INCLUDE_OUTLIERS", true);
  hull_options.flag("CHULL::SEE_NEGLECT", false);
  chull::ConvexHull AEhull;
  AEhull = chull::ConvexHull(hull_options, points);

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " resulting hull has " << AEhull.m_facets.size() << " raw facets" << endl;
  }

  vector<vector<uint>> facet_collection;
  AEhull.getJoinedFacets(facet_collection);
  for (std::vector<vector<uint>>::const_iterator f = facet_collection.begin(); f != facet_collection.end(); ++f) {
    vector<uint> nf;
    for (std::vector<uint>::const_iterator v = f->begin(); v != f->end(); ++v) {
      nf.push_back(*v);
    }
    facets.push_back(nf);
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " after joining " << facets.size() << " facets are remaining" << endl;
  }

  for (std::vector<vector<uint>>::const_iterator f = facets.begin(); f != facets.end(); ++f) {
    if (f->size() < 10) {
      facet_order[f->size() - 3]++;
    } else {
      facet_order[7]++;
    }
  }
  for (size_t t = 0; t < facet_collection.size(); t++) {
    vector<xvector<double>> facet_coords;
    for (std::vector<uint>::const_iterator ind = facet_collection[t].begin(); ind != facet_collection[t].end(); ind++) {
      facet_coords.push_back(points[*ind]);
    }
    facet_area.push_back(aurostd::areaPointsOnPlane(facet_coords));
  }
  volume = aurostd::volume(points, facets, true);
  area = aurostd::sum(facet_area);
  has_hull = true;
}

// ***************************************************************************
// AtomEnvironment::index2Point() - HE20210408
// ***************************************************************************

/// @brief lookup function to map flat neighbor index back into element sorted coordinates_neighbor list
/// @param index neighbor index
/// @return neighbor coordinates
xvector<double> AtomEnvironment::index2Point(uint index) {
  for (size_t i = 0; i < coordinates_neighbor.size(); i++) {
    if (index < coordinations_neighbor[i]) {
      return coordinates_neighbor[i][index];
    } else {
      index -= coordinations_neighbor[i];
    }
  }
  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "index out of bounds", _INDEX_BOUNDS_);
}

// ***************************************************************************
// AtomEnvironment::toJSON() - HE20210408
// ***************************************************************************

/// @brief serialize AtomEnvironment class to json
/// @return json string
/// @authors
/// @mod{HE,20210721,created}
/// @mod{HE,20231029,update to aurostd::JSON}
aurostd::JSON::object AtomEnvironment::toJSON(bool full) const {
  aurostd::JSON::object ae_json(aurostd::JSON::object_types::DICTIONARY);

  // this is the XtalFinder printing mode, we need to retain
  // this printing mode since the work has been published
  // this is flatter than the other printing method //DX20210624
  if (mode == ATOM_ENVIRONMENT_MODE_1) {
    // element_center
    ae_json["element_center"] = element_center;
    // type_center
    ae_json["type_center"] = type_center;
    // elements_neighbor
    ae_json["elements_neighbor"] = elements_neighbor;
    // distances_neighbor
    ae_json["distances_neighbor"] = distances_neighbor;
    // coordinations_neighbor
    ae_json["coordinations_neighbor"] = coordinations_neighbor;
    // coordinates_neighbor
    ae_json["coordinates_neighbor"] = coordinates_neighbor;
  }

  else {
    ae_json["ae_mode"] = mode;
    ae_json["center_element"] = element_center;
    ae_json["center_element_index"] = type_center;
    ae_json["element_count"] = num_types;

    //    "elements_neighbor"
    if (full) {
      ae_json["neighbor_elements"] = aurostd::JSON::object(aurostd::JSON::object_types::LIST);
      for (uint i = 0; i < elements_neighbor.size(); i++) {
        const aurostd::JSON::object distance_element(aurostd::JSON::object_types::DICTIONARY);
        distance_element["index"] = i;
        distance_element["name"] = elements_neighbor[i];
        distance_element["min_distance"] = distances_neighbor[i];
        distance_element["coordination"] = coordinations_neighbor[i];
        ae_json["neighbor_elements"].push_back(distance_element);
      }
    }
    if (has_hull) {
      ae_json["volume"] = volume;
      ae_json["area"] = area;
    }

    ae_json["neighbors"] = aurostd::JSON::object(aurostd::JSON::object_types::LIST);
    uint index = 0;
    for (size_t i = 0; i < coordinations_neighbor.size(); i++) {
      for (uint k = 0; k < coordinations_neighbor[i]; k++) {
        const aurostd::JSON::object neighbor(aurostd::JSON::object_types::DICTIONARY);
        if (full) {
          neighbor["index"] = index;
          neighbor["element"] = elements_neighbor[i];
          neighbor["element_index"] = types_neighbor[i];
          neighbor["coordinate"] = coordinates_neighbor[i][k];
        } else {
          neighbor["element"] = elements_neighbor[i];
          neighbor["coordinate"] = coordinates_neighbor[i][k];
        }
        ae_json["neighbors"].push_back(neighbor);
        index++;
      }
    }

    if (has_hull && full) {
      ae_json["facets"] = aurostd::JSON::object(aurostd::JSON::object_types::LIST);
      for (uint i = 0; i < facets.size(); i++) {
        const aurostd::JSON::object facet_entry(aurostd::JSON::object_types::DICTIONARY);
        facet_entry["area"] = facet_area[i];
        facet_entry["vertices"] = facets[i];
        ae_json["facets"].push_back(facet_entry);
      }

      ae_json["facet_order"] = facet_order;
    }
  }
  return ae_json;
}

// ***************************************************************************
// AtomEnvironment::getAtomEnvironment() - DX20191122
// ***************************************************************************
// determines the atomic environment around a central atom
// current functionality:
//   - calculates the nearest neighbors by element-type, i.e., minimum coordination shell for a given element-type
//     the neighbor elements, types, distance, coordination, and coordinates are stored in the object
//     only one distance is sto
// preliminary functionality, can/will be expanded in the future
void AtomEnvironment::getAtomEnvironment(const xstructure& xstr, uint center_index, uint ae_mode) {
  const vector<string> neighbor_elements;
  getAtomEnvironment(xstr, center_index, neighbor_elements, ae_mode);
}

void AtomEnvironment::getAtomEnvironment(const xstructure& xstr, uint center_index, const vector<string>& neighbor_elements, uint ae_mode) {
  // ---------------------------------------------------------------------------
  // ATOM_ENVIRONMENT_MODE_1 : default minimum coordination shell
  // [FUTURE] ATOM_ENVIRONMENT_MODE_2 : out to a certain radius
  // [FUTURE] ATOM_ENVIRONMENT_MODE_3 : largest gap in radial distribution function

  // ---------------------------------------------------------------------------
  // get central atom info
  mode = ae_mode;
  for (size_t i = 0; i < xstr.atoms.size(); i++) {
    if (i == center_index) {
      element_center = xstr.atoms[i].name;
      type_center = xstr.atoms[i].type;
    }
  }

  num_types = xstr.species.size();

  // ---------------------------------------------------------------------------
  // ATOM_ENVIRONMENT_MODE_1 : minimum coordination environment for each type
  if (mode == ATOM_ENVIRONMENT_MODE_1) {
    for (size_t i = 0; i < xstr.species.size(); i++) {
      // check if types are restricted, otherwise get closest neighbors by type
      if (aurostd::WithinList(neighbor_elements, xstr.species[i]) || neighbor_elements.empty()) {
        uint frequency = 0;
        double min_dist = AUROSTD_MAX_DOUBLE;
        vector<xvector<double>> coordinates;

        // calculate minimum coordination shell to a particular element-type
        minimumCoordinationShell(xstr, center_index, min_dist, frequency, coordinates, xstr.species[i]);

        // store
        elements_neighbor.push_back(xstr.species[i]);
        types_neighbor.push_back(i);
        distances_neighbor.push_back(min_dist);
        coordinations_neighbor.push_back(frequency);
        coordinates_neighbor.push_back(coordinates);
        num_neighbors += frequency;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // [FUTURE] ATOM_ENVIRONMENT_MODE_2 : environment out to a given radius

  // ---------------------------------------------------------------------------
  // [FUTURE] ATOM_ENVIRONMENT_MODE_3 : environment out to largest gap in radial distribution function
  // i.e., GFA convention
}

// ***************************************************************************
// getAtomEnvironments() - DX20191122
// ***************************************************************************
vector<AtomEnvironment> getAtomEnvironments(const xstructure& xstr, uint mode) {
  // Calculate the atomic environments in the structure

  vector<AtomEnvironment> environments;

  for (size_t i = 0; i < xstr.atoms.size(); i++) {
    AtomEnvironment env;
    env.getAtomEnvironment(xstr, i, mode);
    environments.push_back(env);
  }
  return environments;
}

// ***************************************************************************
// writeAtomEnvironments() - HE20210723
// ***************************************************************************

void writeAtomEnvironments(vector<AtomEnvironment> AE, const std::map<string, string> meta_data) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  const aurostd::JSON::object ae_json(aurostd::JSON::object_types::DICTIONARY);
  string file_name = "atomic_environment.json";
  string directory_name;
  string file_path;
  string file_extension = ".json";

  // set filetype
  filetype ftype = json_ft;
  if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
    ftype = json_ft;
    file_extension = ".json";
  } else if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
    ftype = txt_ft;
    file_extension = ".txt";
  }

  // for now just JSON is supported
  if (ftype != json_ft) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Just JSON is supported at the moment.", _INPUT_ERROR_);
  }

  // construct file path
  if (XHOST.vflag_control.flag("FILE")) {
    file_name = XHOST.vflag_control.getattachedscheme("FILE");
  } else {
    if (meta_data.find("auid") != meta_data.end()) {
      const string auid = meta_data.at("auid");
      if (auid.find("aflow:") != std::string::npos) {
        file_name = auid.substr(6);
      } else {
        file_name = auid;
      }
    }
  }

  // ensure that filename has the appropriate extension
  if (!(file_name.size() >= file_extension.size() && 0 == file_name.compare(file_name.size() - file_extension.size(), file_extension.size(), file_extension))) {
    file_name += file_extension;
  };

  if (XHOST.vflag_control.flag("DIRECTORY")) {
    directory_name = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    aurostd::DirectoryMake(directory_name);
  }

  if (!directory_name.empty()) {
    file_path = directory_name + "/" + file_name;
  } else {
    file_path = file_name;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Saving " << AE.size() << " atomic environments" << endl;
  }

  if (!meta_data.empty()) {
    for (const auto& meta_entry : meta_data) {
      ae_json[meta_entry.first] = meta_entry.second;
    }
  }
  ae_json["atomic_environments"] = aurostd::JSON::object(aurostd::JSON::object_types::DICTIONARY);
  for (auto& i : AE) {
    ae_json["atomic_environments"].push_back(i.toJSON());
  }

  aurostd::string2file(ae_json.toString(), file_path, aurostd::compression_type::None, "WRITE");
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Written to " << file_path << endl;
  }
}

// ***************************************************************************
// getLFAAtomEnvironments() - DX20191122
// ***************************************************************************
vector<AtomEnvironment> getLFAAtomEnvironments(const xstructure& xstr, const string& lfa, const vector<string>& LFAs, uint mode) {
  // Calculate the LFA atomic environments
  // i.e., only environments comprised of the least frequent atom (LFA) type
  // use-case: quickly screen for potential isoconfigurational structures (see AFLOW-XtalFinder)

  vector<AtomEnvironment> environments_LFA;

  for (size_t i = 0; i < xstr.atoms.size(); i++) {
    if (xstr.atoms[i].name == lfa) {
      AtomEnvironment LFA_env;
      LFA_env.getAtomEnvironment(xstr, i, LFAs, mode);
      environments_LFA.push_back(LFA_env);
    }
  }
  return environments_LFA;
}

bool sortAtomsTypes(const _atom& a1, const _atom& a2) { // CO20180705
  // CO20190218
  // sorting on DOUBLES is dangerous, we need to avoid flipping equivalent atoms
  // so we need to set a cutoff
  if (a1.type != a2.type) {
    return a1.type < a2.type;
  }
  if (!aurostd::isequal(a1.partial_occupation_value, a2.partial_occupation_value, _AFLOW_POCC_ZERO_TOL_)) {
    return a1.partial_occupation_value > a2.partial_occupation_value;
  } // reverse, we actually want lowest type but highest occupation first
  return a1.basis < a2.basis; // maintain previous relative ordering
  // if(a1.type!=a2.type){return a1.type<a2.type;}
  // double dist1=aurostd::modulus(a1.fpos);
  // double dist2=aurostd::modulus(a2.fpos);
  // return dist1<dist2;
  // prettier (fpos standard)
  // if(a1.fpos.rows!=3 || a1.fpos.rows!=a2.fpos.rows){
  //   cerr << "XSTRUCTURE::sortAtomsNames:: bad cartesian coordinates" << endl;
  //   throw aurostd::xerror(__AFLOW_FILE__,XPID+"sortAtomsTypes():","Throw for debugging purposes.",_GENERIC_ERROR_);
  // }
  // for(uint i=1;i<=3;i++){
  //   if(a1.fpos[i]!=a2.fpos[i]){return a1.fpos[i]<a2.fpos[i];}
  // }
  // return false;
}

// ideal for AFLOW (alphabetized) + POCC (grouped by occupations)
bool sortAtomsNames(const _atom& a1, const _atom& a2) { // CO20180705
  // CO20190218
  // sorting on DOUBLES is dangerous, we need to avoid flipping equivalent atoms
  // so we need to set a cutoff
  if (a1.name != a2.name) {
    return a1.name < a2.name;
  }
  if (!aurostd::isequal(a1.partial_occupation_value, a2.partial_occupation_value, _AFLOW_POCC_ZERO_TOL_)) {
    return a1.partial_occupation_value > a2.partial_occupation_value;
  } // reverse, we actually want lowest type but highest occupation first
  return a1.basis < a2.basis; // maintain previous relative ordering
  // return true;
  // if(a1.name!=a2.name){return a1.name<a2.name;}
  // double dist1=aurostd::modulus(a1.fpos);
  // double dist2=aurostd::modulus(a2.fpos);
  // return dist1<dist2;
  // prettier (fpos standard)
  // if(a1.fpos.rows!=3 || a1.fpos.rows!=a2.fpos.rows){
  //   cerr << "XSTRUCTURE::sortAtomsNames:: bad cartesian coordinates" << endl;
  //   throw aurostd::xerror(__AFLOW_FILE__,XPID+"sortAtomsNames():","Throw for debugging purposes.",_GENERIC_ERROR_);
  // }
  // for(uint i=1;i<=3;i++){
  //   if(a1.fpos[i]!=a2.fpos[i]){return a1.fpos[i]<a2.fpos[i];}
  // }
  // return false;
}

bool sortAtomsDist(const _atom& a1, const _atom& a2) {
  // CO20190218
  // sorting on DOUBLES is dangerous, we need to avoid flipping equivalent atoms
  // so we need to set a cutoff
  if (a1.type != a2.type) {
    return a1.type < a2.type;
  }
  const double dist1 = aurostd::modulus(a1.cpos);
  const double dist2 = aurostd::modulus(a2.cpos);
  if (!aurostd::isequal(dist1, dist2, _ZERO_TOL_)) {
    return dist1 < dist2;
  }
  // CO20180705 - maybe we need to consider adding tol here
  return sortAtomsTypes(a1, a2); // CO20180705, pocc values!
}

// void sortAtomsDist() {std::stable_sort(atoms.begin(),atoms.end(),sortAtomsDist);}

bool sortAtomsEquiv(const _atom& a1, const _atom& a2) {
  if (a1.type != a2.type) {
    return a1.type < a2.type;
  }
  // this is generally implied by equivalent, but not so for POCC, so keep
  if (a1.equivalent != a2.equivalent) {
    return a1.equivalent < a2.equivalent;
  }
  return sortAtomsTypes(a1, a2); // CO 180705, pocc values!
} // CO190101

// ---------------------------------------------------------------------------
// Wyckoff sorting function (by Wyckoff letter) //DX20200515
bool sortWyckoffByLetter(const wyckoffsite_ITC& a, const wyckoffsite_ITC& b) {
  // compare letter
  if (a.letter < b.letter) {
    return true;
  }
  // if letters the same, sort by type
  else if (a.letter == b.letter) {
    if (a.type < b.type) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

// ---------------------------------------------------------------------------
// Wyckoff sorting function (by atom type) //DX20200515
bool sortWyckoffByType(const wyckoffsite_ITC& a, const wyckoffsite_ITC& b) {
  // compare type
  if (a.type < b.type) {
    return true;
  }
  // if types the same, sort by letter
  else if (a.type == b.type) {
    if (a.letter < b.letter) {
      return true;
    } else if (a.letter == b.letter) { // DX20230129 - if the letter is the same
      if (a.parameter_index < b.parameter_index) {
        return true;
      } // DX20230130 - check the parameter index
      else {
        return false;
      } // DX20230130 - otherwise do not change order
    } else {
      return false;
    }
  }
  return false;
}

// DX20180124 - added ibrav to lattice - START
//  **************************************************************************
//  pflow::QE_ibrav2lattice
//  **************************************************************************
namespace pflow {
  xmatrix<double> QE_ibrav2lattice(const int& ibrav, const xvector<double>& parameters, const bool& isabc) {
    // Lattice based on ibrav and lattice parameters (see http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#ibrav)
    const xmatrix<double> lattice;

    double a;
    double b;
    double c;
    double gamma;
    double beta;
    double alpha = 0;

    if (isabc) {
      a = parameters(1); // parameters(1) = a
      b = parameters(2); // parameters(2) = b
      c = parameters(3); // parameters(3) = c
    } else {
      a = parameters(1) * bohr2angstrom; // parameters(1) = a   //DX20180215 - added bohr2angstrom (celldm is in Bohr)
      b = parameters(2) * parameters(1) * bohr2angstrom;
      // parameters(2) = b/a //DX20180215 - added bohr2angstrom (celldm is in Bohr)
      c = parameters(3) * parameters(1) * bohr2angstrom;
      // parameters(3) = c/a //DX20180215 - added bohr2angstrom (celldm is in Bohr)
    }
    gamma = acos(parameters(4));
    // parameters(4) = cos(gamma) NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma
    beta = acos(parameters(5));
    // parameters(5) = cos(beta)  NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma
    alpha = acos(parameters(6));
    // parameters(6) = cos(alpha) NOTE: Different than AFLOW's order convention of a,b,c,alpha,beta,gamma

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

    if (ibrav == 1) { // CUB
      a1 = a * xn;
      a2 = a * yn;
      a3 = a * zn;
    } else if (ibrav == 2) { // FCC
      a1 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * zn;
      a2 = (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
      a3 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn;
    } else if (ibrav == 3) { // BCC
      a1 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
      a2 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
      a3 = -(1.0 / 2.0) * a * xn - (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
    } else if (ibrav == -3) { // BCC (symmetric axis)
      a1 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
      a2 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
      a3 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn - (1.0 / 2.0) * a * zn;
    } else if (ibrav == 4) { // Hexagonal/Trigonal
      a1 = a * xn;
      a2 = -(1.0 / 2.0) * a * xn + (sqrt(3.0) / 2.0) * a * yn;
      a3 = c * zn;
    } else if (ibrav == 5) { // Trigonal R
      a1 = sqrt((1.0 - cos(gamma)) / 2.0) * a * xn - sqrt((1.0 - cos(gamma)) / 6.0) * a * yn + sqrt((1.0 + (2.0 * cos(gamma))) / 3.0) * a * zn;
      a2 = 2.0 * sqrt((1.0 - cos(gamma)) / 6.0) * a * yn + sqrt((1.0 + (2.0 * cos(gamma))) / 3.0) * a * zn;
      a3 = -sqrt((1.0 - cos(gamma)) / 2.0) * a * xn - sqrt((1.0 - cos(gamma)) / 6.0) * a * yn + sqrt((1.0 + (2.0 * cos(gamma))) / 3.0) * a * zn;
    } else if (ibrav == -5) { // Trigonal R (3-fold axis)
      a1 = (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) - 2.0 * sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * xn +
           (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) + sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * yn +
           (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) + sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * zn;
      a2 = (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) + sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * xn +
           (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) - 2.0 * sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * yn +
           (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) + sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * zn;
      a3 = (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) + sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * xn +
           (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) + sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * yn +
           (a / sqrt(3.0)) * ((sqrt((1.0 + (2.0 * cos(gamma))) / 3.0)) - 2.0 * sqrt(2.0) * (sqrt((1.0 - cos(gamma)) / 6.0))) * zn;
    } else if (ibrav == 6) { // TET
      a1 = a * xn;
      a2 = a * yn;
      a3 = c * zn;
    } else if (ibrav == 7) { // BCT
      a1 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * a * yn + (1.0 / 2.0) * c * zn;
      a2 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn + (1.0 / 2.0) * c * zn;
      a3 = -(1.0 / 2.0) * a * xn - (1.0 / 2.0) * a * yn + (1.0 / 2.0) * c * zn;
    } else if (ibrav == 8) { // ORC
      a1 = a * xn;
      a2 = b * yn;
      a3 = c * zn;
    } else if (ibrav == 9) { // ORCC
      a1 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn;
      a2 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn;
      a3 = c * zn;
    } else if (ibrav == -9) { // ORCC - 2
      a1 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * b * yn;
      a2 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn;
      a3 = c * zn;
    } else if (ibrav == 10) { // ORCF
      a1 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * c * zn;
      a2 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn;
      a3 = (1.0 / 2.0) * b * yn + (1.0 / 2.0) * c * zn;
    } else if (ibrav == 11) { // ORCI
      a1 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn + (1.0 / 2.0) * c * zn;
      a2 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn + (1.0 / 2.0) * c * zn;
      a3 = -(1.0 / 2.0) * a * xn - (1.0 / 2.0) * b * yn + (1.0 / 2.0) * c * zn;
    } else if (ibrav == 12) { // MCL (unique axis c)
      a1 = a * xn;
      a2 = cos(gamma) * b * xn + sin(gamma) * b * yn;
      a3 = c * zn;
    } else if (ibrav == -12) { // MCL (unique axis b)
      a1 = a * xn;
      a2 = b * yn;
      a3 = cos(beta) * c * xn + sin(beta) * c * zn;
    } else if (ibrav == 13) { // MCLC
      a1 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * c * zn;
      a2 = cos(gamma) * b * xn + sin(gamma) * b * yn;
      a3 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * c * zn;
    } else if (ibrav == 14) { // TRI
      const double cx = c * cos(beta);
      const double cy = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
      const double cz = sqrt(pow(c, 2.0) - pow(cx, 2.0) - pow(cy, 2.0));

      a1 = a * xn;
      a2 = b * cos(gamma) * xn + b * sin(gamma) * yn;
      a3 = cx * xn + cy * yn + cz * zn;
    }

    lattice(1, 1) = a1(1);
    lattice(1, 2) = a1(2);
    lattice(1, 3) = a1(3);
    lattice(2, 1) = a2(1);
    lattice(2, 2) = a2(2);
    lattice(2, 3) = a2(3);
    lattice(3, 1) = a3(1);
    lattice(3, 2) = a3(2);
    lattice(3, 3) = a3(3);

    return lattice;
  }
} // namespace pflow

// DX20180124 - added ibrav to lattice - START

// **************************************************************************
// Function GetVols, det,
// **************************************************************************
// Wrap up for GetVol functions ... useful to bring
// convasp framework to aflow.

double GetVol(const xmatrix<double>& lat) {
  return std::abs(det(lat));
} // AFLOW_FUNCTION_IMPLEMENTATION
double det(const xvector<double>& v1, const xvector<double>& v2, const xvector<double>& v3) {
  return scalar_product(v1, vector_product(v2, v3));
}

double det(const double& a11, const double& a12, const double& a13, const double& a21, const double& a22, const double& a23, const double& a31, const double& a32, const double& a33) {
  return a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31 - a12 * a21 * a33 - a11 * a23 * a32;
}

double GetVol(const xvector<double>& v1, const xvector<double>& v2, const xvector<double>& v3) {
  return std::abs(det(v1, v2, v3));
} // AFLOW_FUNCTION_IMPLEMENTATION

// ***************************************************************************
// Function Getabc_angles
// ***************************************************************************
// This function returns a,b,c,alpha,beta,gamma for the
// cell given the lattice xvectors. Dane Morgan, adjusted by SC

// #define _Getabc_angles Getabc_angles
// #define _Getabc_angles __NO_USE_Sortabc_angles

xvector<double> Getabc_angles(const xmatrix<double>& lat, int mode) { // AFLOW_FUNCTION_IMPLEMENTATION
  const xvector<double> data(6);
  data(1) = aurostd::modulus(lat(1));
  data(2) = aurostd::modulus(lat(2));
  data(3) = aurostd::modulus(lat(3));
  data(4) = angle(lat(2), lat(3));
  data(5) = angle(lat(3), lat(1));
  data(6) = angle(lat(1), lat(2));
  if (mode == DEGREES) {
    data(4) *= rad2deg;
    data(5) *= rad2deg;
    data(6) *= rad2deg;
  }
  return data;
}

xvector<double> Getabc_angles(const xmatrix<double>& lat, const xvector<int>& permut, int mode) {
  // AFLOW_FUNCTION_IMPLEMENTATION
  const xvector<double> data(6);
  data(1) = aurostd::modulus(lat(1));
  data(2) = aurostd::modulus(lat(2));
  data(3) = aurostd::modulus(lat(3));
  data(4) = angle(lat(2), lat(3));
  data(5) = angle(lat(3), lat(1));
  data(6) = angle(lat(1), lat(2));
  if (mode == DEGREES) {
    data(4) *= rad2deg;
    data(5) *= rad2deg;
    data(6) *= rad2deg;
  }
  // with permutation - from AVDV
  int i;
  int imin;
  int imax;
  int imid;
  double dmin;
  double dmax;
  dmax = dmin = data(1);
  imax = imin = imid = 1;
  for (i = 1; i <= 3; i++) {
    if (dmax <= data(i)) {
      dmax = data(i);
      imax = i;
    }
    if (dmin > data(i)) {
      dmin = data(i);
      imin = i;
    }
  }
  // set imid
  for (i = 1; i <= 3; i++) {
    if (i != imin && i != imax) {
      imid = i;
    }
  }
  // if all lattice parameters are equal length, numerical noise may cause imin=imax
  if (imin == imax) {
    for (i = 1; i <= 3; i++) {
      if (i != imin && i != imid) {
        imax = i;
      }
    }
  }
  permut[1] = imin;
  permut[2] = imid;
  permut[3] = imax;
  // done permutation - from AVDV
  return data;
}

xvector<double> Getabc_angles(const xvector<double>& r1, // AFLOW_FUNCTION_IMPLEMENTATION
                              const xvector<double>& r2, // AFLOW_FUNCTION_IMPLEMENTATION
                              const xvector<double>& r3, // AFLOW_FUNCTION_IMPLEMENTATION
                              int mode) { // AFLOW_FUNCTION_IMPLEMENTATION
  const xmatrix<double> lat(3, 3);
  lat(1, 1) = r1(1);
  lat(1, 2) = r1(2);
  lat(1, 3) = r1(3);
  lat(2, 1) = r2(1);
  lat(2, 2) = r2(2);
  lat(2, 3) = r2(3);
  lat(3, 1) = r3(1);
  lat(3, 2) = r3(2);
  lat(3, 3) = r3(3);
  return Getabc_angles(lat, mode);
}

xvector<double> Getabc_angles(const xvector<double>& r1, // AFLOW_FUNCTION_IMPLEMENTATION
                              const xvector<double>& r2, // AFLOW_FUNCTION_IMPLEMENTATION
                              const xvector<double>& r3, // AFLOW_FUNCTION_IMPLEMENTATION
                              const xvector<int>& permut, // AFLOW_FUNCTION_IMPLEMENTATION
                              int mode) { // AFLOW_FUNCTION_IMPLEMENTATION
  const xmatrix<double> lat(3, 3);
  lat(1, 1) = r1(1);
  lat(1, 2) = r1(2);
  lat(1, 3) = r1(3);
  lat(2, 1) = r2(1);
  lat(2, 2) = r2(2);
  lat(2, 3) = r2(3);
  lat(3, 1) = r3(1);
  lat(3, 2) = r3(2);
  lat(3, 3) = r3(3);
  return Getabc_angles(lat, permut, mode);
}

xvector<double> Sortabc_angles(const xmatrix<double>& lat, int mode) { // AFLOW_FUNCTION_IMPLEMENTATION
  // with permutation - from AVDV
  int i;
  int imin;
  int imax;
  int imid;
  double dmin;
  double dmax;
  dmax = dmin = aurostd::modulus(lat(1));
  imax = imin = imid = 1;
  for (i = 1; i <= 3; i++) {
    if (dmax <= aurostd::modulus(lat(i))) {
      dmax = aurostd::modulus(lat(i));
      imax = i;
    }
    if (dmin > aurostd::modulus(lat(i))) {
      dmin = aurostd::modulus(lat(i));
      imin = i;
    }
  }
  // set imid
  for (i = 1; i <= 3; i++) {
    if (i != imin && i != imax) {
      imid = i;
    }
  }
  // if all lattice parameters are equal length, numerical noise may cause imin=imax
  if (imin == imax) {
    for (i = 1; i <= 3; i++) {
      if (i != imin && i != imid) {
        imax = i;
      }
    }
  }

  const xvector<double> data(6);
  data(1) = aurostd::modulus(lat(imin));
  data(2) = aurostd::modulus(lat(imid));
  data(3) = aurostd::modulus(lat(imax));
  data(4) = angle(lat(imid), lat(imax));
  data(5) = angle(lat(imin), lat(imax));
  data(6) = angle(lat(imin), lat(imid));
  if (mode == DEGREES) {
    data(4) *= rad2deg;
    data(5) *= rad2deg;
    data(6) *= rad2deg;
  }
  return data;
}

// **************************************************************************
// Function GetClat
// **************************************************************************
// This function returns cartesian lattice xvectors for the cell given
// a,b,c,alpha,beta,gamma. Assumes angles are in radians.
// This routine aligns a along +X, b in XY plane with
// angle gamma from a in +Y direction (bx=b.a_unit=b*cos[gamma],
// by=sqrt(b^2-bx^2)=b*sin[gamma]), and c is then given by
// cx=c.a_unit=c*cos[beta],
// cy=c.by_unit=c*(cos[alpha]-cos[gamma]*cos[beta])/sin[gamma],
// cz=sqrt(c^2-cx^2-cy^2)
// Dane Morgan, adjusted by SC

xmatrix<double> GetClat(const xvector<double>& abc_angles) { // AFLOW_FUNCTION_IMPLEMENTATION
  stringstream message;
  const xmatrix<double> clattice(3, 3);
  const double a = abc_angles[1];
  const double b = abc_angles[2];
  const double c = abc_angles[3];
  const double bc = abc_angles[4] * deg2rad; // angle from b to c (remove a)
  const double ca = abc_angles[5] * deg2rad; // angle from c to a (remove b)
  const double ab = abc_angles[6] * deg2rad; // angle from a to b (remove c)
  //  if(abs(bc)>6.3||abs(ca)>6.3||abs(ab)>6.3) { cerr << _AUROSTD_XLIBS_ERROR_ << "GetClat: angles must be in RADIANT " << endl;throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);}
  clattice(1, 1) = a;
  clattice(2, 1) = b * cos(ab);
  clattice(2, 2) = b * sin(ab);
  clattice(3, 1) = c * cos(ca);
  if (ab < 0.00000001) {
    message << "The angle gamma from a to b is too small" << endl;
    message << "gamma = " << ab << endl;
    message << "STOPPING " << endl;
    message << _AUROSTD_XLIBS_ERROR_ << "ERROR: STOPPING " << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
  }
  clattice(3, 2) = c * (cos(bc) - cos(ab) * cos(ca)) / sin(ab);
  clattice(3, 3) = sqrt(std::abs(c * c - clattice(3, 2) * clattice(3, 2) - clattice(3, 1) * clattice(3, 1)));
  return clattice;
}

xmatrix<double> GetClat(const double& a, const double& b, const double& c, const double& alpha, const double& beta, const double& gamma) {
  stringstream message;
  const xmatrix<double> clattice(3, 3);
  const double bc = alpha * deg2rad; // angle from b to c (remove a)
  const double ca = beta * deg2rad; // angle from c to a (remove b)
  const double ab = gamma * deg2rad; // angle from a to b (remove c)
  //  if(abs(bc)>6.3||abs(ca)>6.3||abs(ab)>6.3) { cerr << _AUROSTD_XLIBS_ERROR_ << "GetClat: angles must be in RADIANT " << endl;throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);}
  clattice(1, 1) = a;
  clattice(2, 1) = b * cos(ab);
  clattice(2, 2) = b * sin(ab);
  clattice(3, 1) = c * cos(ca);
  if (ab < 0.00000001) {
    message << "The angle gamma from a to b is too small" << endl;
    message << "gamma = " << ab << endl;
    message << "STOPPING " << endl;
    message << _AUROSTD_XLIBS_ERROR_ << "ERROR: STOPPING " << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
  }
  clattice(3, 2) = c * (cos(bc) - cos(ab) * cos(ca)) / sin(ab);
  clattice(3, 3) = sqrt(std::abs(c * c - clattice(3, 2) * clattice(3, 2) - clattice(3, 1) * clattice(3, 1)));
  return clattice;
}

// **************************************************************************
// Function GetIntpolStr
// **************************************************************************
// This function gets the structure that is linearly interpolated
// a fraction f of the way between strA and strB.  Interpolation
// nincludes lattice parameters.  The scale factors of the interpolated
// structures are all set to 1.

xstructure GetIntpolStr(xstructure strA, xstructure strB, const double& f, const string& path_flag) {
  strA = ReScale(strA, 1.0);
  strB = ReScale(strB, 1.0);
  // Get new lattice params.
  const xmatrix<double> lati(3, 3);
  xmatrix<double> latA(3, 3);
  latA = strA.lattice;
  xmatrix<double> latB(3, 3);
  latB = strB.lattice;
  xvector<double> dl(3);
  for (int ic = 1; ic <= 3; ic++) {
    dl = f * (latB(ic) - latA(ic));
    for (int jc = 1; jc <= 3; jc++) {
      lati(ic, jc) = latA(ic, jc) + dl(jc);
    }
  }
  // Get new cart. coords.
  if (strA.atoms.size() != strB.atoms.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, _AUROSTD_XLIBS_ERROR_ + " number of atoms must be the same in both structures!!", _INPUT_ILLEGAL_);
  }
  const int size = strA.atoms.size();
  vector<xvector<double>> cposi(size, 3);
  for (int iat = 0; iat < size; iat++) {
    xvector<double> dp(3);
    dp = strB.atoms.at(iat).cpos - strA.atoms.at(iat).cpos;
    // If path_flag is n/N then the path is taken between
    // nearest images.  Otherwise path is between the atoms given.
    if (path_flag == "n" || path_flag == "N") {
      xvector<double> ddp(3);
      ddp = C2F(lati, dp);
      for (int ic = 1; ic <= 3; ic++) {
        ddp(ic) = ddp(ic) - nint(ddp(ic));
      }
      dp = F2C(lati, ddp);
    }
    dp = f * dp;
    cposi.at(iat) = strA.atoms.at(iat).cpos + dp;
  }
  xstructure stri = strA;
  stri.lattice = lati;
  for (int iat = 0; iat < size; iat++) {
    stri.atoms.at(iat).cpos = cposi.at(iat);
    stri.atoms.at(iat).fpos = C2F(stri.lattice, stri.atoms.at(iat).cpos);
  }
  return stri;
}

// **************************************************************************
// Function RadiusSphereLattice
// **************************************************************************
// This function returns the radius of the sphere encompassing the whole cell
double RadiusSphereLattice(const xmatrix<double>& lattice, double scale) {
  double radius = 0;
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        if (aurostd::modulus((double) i * lattice(1) + (double) j * lattice(2) + (double) k * lattice(3)) > radius) {
          radius = aurostd::modulus((double) i * lattice(1) + (double) j * lattice(2) + (double) k * lattice(3));
        }
      }
    }
  }
  return scale * radius;
}

// **************************************************************************
// Function LatticeDimensionSphere
// **************************************************************************
// This function, given a lattice with vectors lat[3][3], finds the
// dimensions along the unit cell vectors such that a sphere of given radius
// fits within a uniform grid of 2dim[1]x2dim[2]x2dim[3] lattice points
// centered at the origin.
//
// The algorithm works by getting the normal (e.g. n1) to each pair of lattice
// vectors (e.g. a2, a3), scaling this normal to have length radius and
// then projecting this normal parallel to the a2,a3 plane onto the
// remaining lattice vector a1. This will tell us the number of a1 vectors
// needed to make a grid to encompass the sphere.
//
// Written as lat_dimension by Anton Van Der Ven
// Adjusted by Stefano Curtarolo

xvector<int> LatticeDimensionSphere(const xmatrix<double>& _lattice, double radius, double scale) {
  // Adapted from AVDV routine
#if DEBUG_LATTICE_DIMENSIONS
  bool LDEBUG = (false || XHOST.DEBUG);
#endif
  xmatrix<double> lattice;
  lattice = scale * _lattice;
  int i;
  int j;
  int k;
  xmatrix<double> invlattice(3, 3);
  xmatrix<double> normals(3, 3);
  const xmatrix<double> frac_normals(3, 3);
  const xvector<int> dim(3);
  double length;
  normals.clear();
  // get the normals to pairs of lattice vectors of length radius
  for (int m = 1; m <= 3; m++) {
    for (int n = 1; n <= 3; n++) {
      for (int l = 1; l <= 3; l++) {
        normals(1, l) += aurostd::eijk(l, m, n) * lattice(2, m) * lattice(3, n);
        normals(2, l) += aurostd::eijk(l, m, n) * lattice(3, m) * lattice(1, n);
        normals(3, l) += aurostd::eijk(l, m, n) * lattice(1, m) * lattice(2, n);
      }
    }
  }
#if DEBUG_LATTICE_DIMENSIONS
  if (LDEBUG) {
    for (uint i = 1; i < (uint) normals.rows + 1; i++) {
      cerr << __AFLOW_FUNC__ << " normals(" << i << ")=" << normals(i) << endl;
    }
  }
#endif
  for (i = 1; i <= 3; i++) {
    length = aurostd::modulus(normals(i));
    for (j = 1; j <= 3; j++) {
      normals(i, j) = radius * normals(i, j) / length;
    }
  }
  // get the normals in the coordinates system of the lattice vectors
  invlattice = aurostd::inverse(lattice);
#if DEBUG_LATTICE_DIMENSIONS
  if (LDEBUG) { // CO20190520
    cerr << __AFLOW_FUNC__ << " normals=" << endl;
    cerr << normals << endl; // CO20190520
    cerr << __AFLOW_FUNC__ << " lattice=" << endl;
    cerr << lattice << endl; // CO20190520
    cerr << __AFLOW_FUNC__ << " inverse(lattice)=" << endl;
    cerr << invlattice << endl; // CO20190520
  } // CO20190520
#endif

  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++) {
      frac_normals(i, j) = 0.0;
      for (k = 1; k <= 3; k++) {
        frac_normals(i, j) = frac_normals(i, j) + normals(i, k) * invlattice(k, j);
      }
    }
  }
  // the diagonals of frac_normal contain the dimensions of the lattice grid that
  // encompasses a sphere of radius = radius
  for (i = 1; i <= 3; i++) {
#if DEBUG_LATTICE_DIMENSIONS
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " abs(frac_normals(i,i))=" << std::abs(frac_normals(i, i)) << endl;
    }
#endif
    dim(i) = (int) ceil(std::abs(frac_normals(i, i)));
  }
  if (max(dim) == 0) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, _AUROSTD_XLIBS_ERROR_ + " dim=0!!", _INPUT_ILLEGAL_);
  }
  return dim;
}

xvector<int> LatticeDimensionSphere(const xstructure& str, double radius) {
  return LatticeDimensionSphere(str.lattice, radius, str.scale); // str.scale*str.lattice,radius);  //CO20171024
}

// **************************************************************************
// F2C and C2F exchange transformations
// **************************************************************************
// Stefano Curtarolo
// xf=(CF)*xc and xc=(FC)*xf   (FC)=inv(CF)
//  A is a property so it is a 2 indices tensor
//  x'=A*x and xf'=Af*xf and xc'=Ac*xc
//  xf'=(CF)*xc'=Af*xf=Af*(CF)*xc   mult inv(CF) *
//  xc'=inv(CF)*Af*(CF)*xc HENCE
//  Ac=inv(CF)*Af*(CF);
//  Af=inv(FC)*Ac*(FC)
//  in another way
//  Ac_ij=inv(CF)_il Af_lm (CF)_mj
//  Af_ij=inv(FC)_il Ac_lm (FC)_mj
//  same for operations, although here FC=inv(CF)
//  see notes on operation of force constant tensors

// **************************************************************************
// -------------------------------------------------------- for COLUMN vectors
xvector<double> F2C(const double& scale, const xmatrix<double>& lattice, const xvector<double>& fpos) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(lattice);
  return f2c * fpos; // fpos are F components per COLUMS !
}

xvector<double> F2C(const xmatrix<double>& lattice, const xvector<double>& fpos) {
  return F2C(1.0, lattice, fpos); // VASP cartesian coordinates are intended normalized on the scale
}

xvector<double> C2F(const double& scale, const xmatrix<double>& lattice, const xvector<double>& cpos) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(lattice);
  return inverse(f2c) * cpos; // cpos are C components per COLUMS !
}

xvector<double> C2F(const xmatrix<double>& lattice, const xvector<double>& cpos) {
  return C2F(1.0, lattice, cpos); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// ---------------------------------------------- for a set of COLUM matrices
xmatrix<double> F2C(const double& scale, const xmatrix<double>& lattice, const xmatrix<double>& fpos) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(lattice);
  return f2c * fpos; // fpos are F components per COLUMS !
}

xmatrix<double> F2C(const xmatrix<double>& lattice, const xmatrix<double>& fpos) {
  return F2C(1.0, lattice, fpos); // VASP cartesian coordinates are intended normalized on the scale
}

xmatrix<double> C2F(const double& scale, const xmatrix<double>& lattice, const xmatrix<double>& cpos) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(lattice);
  return inverse(f2c) * cpos; // cpos are C components per COLUMS !
}

xmatrix<double> C2F(const xmatrix<double>& lattice, const xmatrix<double>& cpos) {
  return C2F(1.0, lattice, cpos); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// -------------------------------------------------------- for _atoms with lattice
_atom F2C(const double& scale, const xmatrix<double>& lattice, const _atom& iatom) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(lattice);
  _atom oatom(iatom);
  oatom.cpos = f2c * iatom.fpos;
  return oatom;
}

_atom F2C(const xmatrix<double>& lattice, const _atom& iatom) {
  return F2C(1.0, lattice, iatom); // VASP cartesian coordinates are intended normalized on the scale
}

_atom C2F(const double& scale, const xmatrix<double>& lattice, const _atom& iatom) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(lattice);
  _atom oatom(iatom);
  oatom.fpos = inverse(f2c) * iatom.cpos;
  return oatom;
}

_atom C2F(const xmatrix<double>& lattice, const _atom& iatom) {
  return C2F(1.0, lattice, iatom); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// -------------------------------------------------------- for _atoms with structure
_atom F2C(const double& scale, const xstructure& str, const _atom& iatom) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(str.lattice);
  _atom oatom(iatom);
  oatom.cpos = f2c * iatom.fpos;
  return oatom;
}

_atom F2C(const xstructure& str, const _atom& iatom) {
  return F2C(1.0, str.lattice, iatom); // VASP cartesian coordinates are intended normalized on the scale
}

_atom C2F(const double& scale, const xstructure& str, const _atom& iatom) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(str.lattice);
  _atom oatom(iatom);
  oatom.fpos = inverse(f2c) * iatom.cpos;
  return oatom;
}

_atom C2F(const xstructure& str, const _atom& iatom) {
  return C2F(1.0, str.lattice, iatom); // VASP cartesian coordinates are intended normalized on the scale
}

// **************************************************************************
// **************************************************************************
// -------------------------------------------------------- for OPERATORS (matrices)
xmatrix<double> FF2CC(const double& scale, const xmatrix<double>& lattice, const xmatrix<double>& fmat) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(lattice);
  return f2c * fmat * inverse(f2c); // fmat is an operation in F coordinates
}

xmatrix<double> FF2CC(const xmatrix<double>& lattice, const xmatrix<double>& fmat) {
  return FF2CC(1.0, lattice, fmat); // VASP cartesian coordinates are intended normalized on the scale
}

xmatrix<double> CC2FF(const double& scale, const xmatrix<double>& lattice, const xmatrix<double>& cmat) {
  xmatrix<double> f2c(3, 3);
  f2c = scale * trasp(lattice);
  return inverse(f2c) * cmat * f2c; // cmat is an operation in C coordinates
}

xmatrix<double> CC2FF(const xmatrix<double>& lattice, const xmatrix<double>& cmat) {
  return CC2FF(1.0, lattice, cmat); // VASP cartesian coordinates are intended normalized on the scale
}

// ***************************************************************************
// Function IdenticalAtoms
// ***************************************************************************
// Makes all atoms the same. Stefano Curtarolo Nov 2008
void xstructure::IdenticalAtoms() {
  xvector<double> fpos(3);
  xvector<double> cpos(3);
  // fix atoms
  for (size_t i = 0; i < atoms.size(); i++) {
    fpos = atoms[i].fpos; // save fpos
    cpos = atoms[i].cpos; // save cpos
    atoms[i] = atoms.at(0); // use this so to copy all info
    atoms[i].fpos = fpos; // fix back fpos
    atoms[i].cpos = cpos; // fix back cpos
    atoms[i].partial_occupation_flag = false;
    atoms[i].partial_occupation_value = 1.0;
  }
  // fix numbers
  num_each_type.clear();
  num_each_type.push_back(atoms.size());
  comp_each_type.clear();
  comp_each_type.push_back((double) atoms.size());
  GetStoich(); // CO20170724
  // fix species
  for (size_t i = 1; i < species.size(); i++) {
    species[i] = "";
    species_pp.at(i) = "";
    species_volume.at(i) = 0.0;
    species_mass.at(i) = 0.0;
  }
}

xstructure IdenticalAtoms(const xstructure& str) {
  xstructure sstr(str);
  sstr.IdenticalAtoms();
  return sstr;
}

// ***************************************************************************
// Function SwapCoordinates
// ***************************************************************************
// Permute Coordinates i with j Stefano Curtarolo Oct 2009
// Wahyu Setyawan: keep the same right-handness of the original
void xstructure::SwapCoordinates(const uint& ii, const uint& jj) { // Permute Coordinates i with j
  // is obvious
  uint kk = 0;
  if ((ii == 1 && jj == 2) || (ii == 2 && jj == 1)) {
    kk = 3; // permutation
  }
  if ((ii == 2 && jj == 3) || (ii == 3 && jj == 2)) {
    kk = 1; // permutation
  }
  if ((ii == 3 && jj == 1) || (ii == 1 && jj == 3)) {
    kk = 2; // permutation
  }
  if ((ii == jj) || (ii < 1) || (ii > 3) || (jj < 1) || (jj > 3)) {
    return; // nothing to do
  }

  // do the lattice
  double tmp;
  for (uint i = 1; i <= 3; i++) {
    tmp = lattice[ii][i];
    lattice[ii][i] = lattice[jj][i];
    lattice[jj][i] = tmp; // swap RAW_ii with RAW_jj
    lattice[kk][i] = -lattice[kk][i]; // keep the right-handness of the original
  }
  FixLattices(); // fix the lattices ...
  // cartesian atoms are the same (crystal was not really changed, only definitions of abc were
  // so now I can get the new fractional
  for (size_t i = 0; i < atoms.size(); i++) {
    atoms[i].fpos = C2F(lattice, atoms[i].cpos); // DONE, now it has new fpos with the same cpos !
  }
  // done bring all in cell now
  BringInCell();
}

xstructure SwapCoordinates(const xstructure& str, const uint& ii, const uint& jj) {
  xstructure sstr(str);
  sstr.SwapCoordinates(ii, jj);
  return sstr;
}

// ***************************************************************************
// Function GetLatticeType()
// ***************************************************************************
// Determine the real, reciprocal, and superlattice symmetry information
// Stefano Curtarolo
// Modified by David Hicks (DX)
// Includes self-consistency loop to ensure descriptions are commensurate
// DX20210225 - cleaned/consolidated function
void xstructure::GetLatticeType(double sym_eps, bool no_scan) {
  xstructure str_sp;
  xstructure str_sc;
  GetLatticeType(str_sp, str_sc, sym_eps, no_scan);
}

void xstructure::GetLatticeType(xstructure& str_sp, xstructure& str_sc, double sym_eps, bool no_scan) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if (tolerance == AUROSTD_MAX_DOUBLE) {
    if ((*this).sym_eps_calculated) {
      tolerance = (*this).sym_eps;
    } else {
      tolerance = SYM::defaultTolerance((*this));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl;
  }

  // keep track of self-consistent tolerance
  bool same_eps = false;
  uint count = 0;
  const uint count_max = 100; // safety for while loop, don't calculate forever

  // update tolerance info in *this
  (*this).sym_eps = tolerance;
  (*this).sym_eps_calculated = true;
  const double tolerance_orig = tolerance; // DX20210623

  // ---------------------------------------------------------------------------
  // loop over the real, reciprocal, and superlattice analysis until all
  // symmetries are commensurate at a common tolerance value
  while (!same_eps && count++ < count_max) {
    // ---------------------------------------------------------------------------
    // clear to start
    (*this).ClearSymmetry();

    // ---------------------------------------------------------------------------
    // update the tolerance, it may have change during loop
    tolerance = (*this).sym_eps;
    no_scan = (*this).sym_eps_no_scan;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [2] Top of self-consistent lattice-type loop (calculating real, reciprocal, and superlattice types) (sym_eps=" << tolerance
           << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl;
    }

    // ---------------------------------------------------------------------------
    // check if consistency checks failed (maxed while loop iteration)
    // turn of scan
    if (count == count_max) {
      no_scan = (*this).sym_eps_no_scan = true;
      tolerance = tolerance_orig;
      // set to original tolerance //DX20210623 - originally sym_eps, but this could be AUROSTD_MAX_DOUBLE;
      cerr << __AFLOW_FUNC__ << " Unable to calculate consistent symmetry. Calculating at original tolerance (sym_eps=" << sym_eps << ") and ignoring consistency checks." << endl;
    }

    // ---------------------------------------------------------------------------
    // REAL - pass in str_sp and str_sc to keep primitive and conventional info
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [3] Calculate real lattice type (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl;
    }
    (*this).GetRealLatticeType(str_sp, str_sc, tolerance);
    tolerance = (*this).sym_eps; // update the tolerance

    // ---------------------------------------------------------------------------
    // RECIPROCAL
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [4] Calculate reciprocal lattice type (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl;
    }
    (*this).GetReciprocalLatticeType(tolerance);
    if (!no_scan && (*this).sym_eps != tolerance) {
      continue;
    } // if tolerance changed, recalc

    // ---------------------------------------------------------------------------
    // SUPERLATTICE
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [5] Calculate superlattice type (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl;
    }
    (*this).GetSuperlatticeType(tolerance);
    if (!no_scan && (*this).sym_eps != tolerance) {
      continue;
    } // if tolerance changed, recalc

    // made it to the end with same sym_eps
    same_eps = true;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [6] Lattice types calculation finished! (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl;
  }
}

// ***************************************************************************
// Function GetExtendedCrystallographicData() //DX20210302
// ***************************************************************************
// Determine the real, reciprocal, superlattice, and space group symmetries
// Includes self-consistency loop to ensure descriptions are commensurate
void xstructure::GetExtendedCrystallographicData(double sym_eps, bool no_scan, int setting) {
  xstructure str_sp;
  xstructure str_sc;
  GetExtendedCrystallographicData(str_sp, str_sc, sym_eps, no_scan, setting);
}

void xstructure::GetExtendedCrystallographicData(xstructure& str_sp, xstructure& str_sc, double sym_eps, bool no_scan, int setting) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if (tolerance == AUROSTD_MAX_DOUBLE) {
    if ((*this).sym_eps_calculated) {
      tolerance = (*this).sym_eps;
    } else {
      tolerance = SYM::defaultTolerance((*this));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl;
  }

  // keep track of self-consistent tolerance
  const bool force_perform = true;
  bool same_eps = false;
  uint count = 0;
  const uint count_max = 100; // safety for while loop, don't calculate forever

  // update tolerance info in *this
  (*this).sym_eps = tolerance;
  (*this).sym_eps_calculated = true;
  const double tolerance_orig = tolerance; // DX20210623

  // ---------------------------------------------------------------------------
  // loop over the real, reciprocal, and superlattice analysis until all
  // symmetries are commensurate with a common tolerance value
  while (!same_eps && count++ < count_max) {
    // ---------------------------------------------------------------------------
    // clear to start
    (*this).ClearSymmetry();

    // ---------------------------------------------------------------------------
    // update the tolerance, it may have change during loop
    tolerance = (*this).sym_eps;
    no_scan = (*this).sym_eps_no_scan;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [2] Top of self-consistent extended crystallographic data loop (calculating lattice type and space group data) (sym_eps=" << tolerance
           << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl;
    }

    // ---------------------------------------------------------------------------
    // check if consistency checks failed (maxed while loop iteration)
    // turn off scan
    if (count == count_max) {
      no_scan = (*this).sym_eps_no_scan = true;
      tolerance = tolerance_orig;
      // set to original tolerance //DX20210623 - originally sym_eps, but this could be AUROSTD_MAX_DOUBLE;
      cerr << __AFLOW_FUNC__ << " Unable to calculate consistent symmetry. Calculating at original tolerance (sym_eps=" << sym_eps << ") and ignoring consistency checks." << endl;
    }

    // ---------------------------------------------------------------------------
    // REAL, RECIPROCAL, and SUPERLATTICE data
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [3] Calculate real, reciprocal, and superlattice information (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl;
    }
    (*this).GetLatticeType(str_sp, str_sc, (*this).sym_eps);
    tolerance = (*this).sym_eps; // update the tolerance

    // ---------------------------------------------------------------------------
    // space group data
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [4] Calculate the space group symmetry information (sym_eps=" << tolerance << ", sym_eps_change_count=" << (*this).sym_eps_change_count << ")" << endl;
    }
    (*this).SpaceGroup_ITC(tolerance, -1, setting, no_scan);
    if (!no_scan && (*this).sym_eps != tolerance) {
      continue;
    } // if tolerance changed, recalc

    // ---------------------------------------------------------------------------
    // check GetLatticeType vs SpaceGroup_ITC results
    if (!(*this).sym_eps_no_scan) {
      int multiplicity_of_primitive = str_sp.fgroup.size() / str_sp.pgroup_xtal.size();
      bool derivative_structure = false;
      const string lattice_and_centering = LATTICE::Lattice2TypeAndCentering((*this).bravais_lattice_type);
      // DX20210412 - check centering
      const string lattice_and_centering_from_sg = SYM::spacegroup2latticeAndCentering((*this).space_group_ITC);
      // DX20210412 - check centering
      if (!(lattice_and_centering == lattice_and_centering_from_sg && SYM::ComparePointGroupAndSpaceGroupString((*this), multiplicity_of_primitive, derivative_structure))) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " WARNING: Space group symbol and point group symbol do not match. (sg=" << GetSpaceGroupName((*this).space_group_ITC, (*this).directory)
               << ", centering_sg=" << lattice_and_centering_from_sg << " | pg=" << (*this).point_group_Hermann_Mauguin << ", centering=" << lattice_and_centering << ") [dir=" << (*this).directory << "]" << endl;
        }
        if (!SYM::change_tolerance((*this), (*this).sym_eps, (*this).dist_nn_min, (*this).sym_eps_no_scan)) {
          if (force_perform) {
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " WARNING: Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << (*this).directory << endl;
            }
            no_scan = (*this).sym_eps_no_scan = true; // DX20210331
            (*this).sym_eps = tolerance_orig;
            // set to original tolerance //DX20210623 - originally sym_eps, but this could be AUROSTD_MAX_DOUBLE;
          }
        }
        continue;
      }
    }

    // made it to the end with same sym_eps
    same_eps = true;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [5] Extended crystallographic data calculation finished! (sym_eps=" << tolerance << ", sym_eps_change_count=" << count << ")" << endl;
  }
}

// ***************************************************************************
// Function GetRealLatticeType //DX20210209
// ***************************************************************************
// Determine the real space symmetry information
// Includes self-consistency loop to ensure descriptions are commensurate
void xstructure::GetRealLatticeType(double sym_eps) {
  xstructure str_sp;
  xstructure str_sc;
  GetRealLatticeType(str_sp, str_sc, sym_eps);
}

void xstructure::GetRealLatticeType(xstructure& str_sp, xstructure& str_sc, double sym_eps) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if (tolerance == AUROSTD_MAX_DOUBLE) {
    if ((*this).sym_eps_calculated) {
      tolerance = (*this).sym_eps;
    } else {
      tolerance = SYM::defaultTolerance((*this));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl;
  }

  // need to create a copy
  xstructure str_in = *this;

  // update tolerance info for all xstructures
  str_in.sym_eps = str_sp.sym_eps = str_sc.sym_eps = tolerance;
  str_in.sym_eps_calculated = str_sp.sym_eps_calculated = str_sc.sym_eps_calculated = (*this).sym_eps_calculated;
  str_in.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = (*this).sym_eps_change_count;
  str_in.sym_eps_no_scan = str_sp.sym_eps_no_scan = str_sc.sym_eps_no_scan = (*this).sym_eps_no_scan; // DX20210430

  // calculate
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [2]" << endl;
  }
  LATTICE::Bravais_Lattice_StructureDefault(str_in, str_sp, str_sc); // STD tolerance  // ONLY BRAVAIS_CRYSTAL

  // set properties
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [3]" << endl;
  }
  if (str_sp.pgroup_calculated == false) {
    str_sp.CalculateSymmetryPointGroup(false); // cerr << "POINT GROUP" << endl;
  }
  if (str_sp.fgroup_calculated == false) {
    str_sp.CalculateSymmetryFactorGroup(false); // cerr << "FACTOR GROUP" << endl;
  }
  if (str_sp.pgroup_xtal_calculated == false) {
    str_sp.CalculateSymmetryPointGroupCrystal(false);
  }
  // cerr << "POINT GROUP XTAL" << endl;
  //   *this=str_sp; // more obvious but will mess up the structures.... we only want to take the properties
  this->bravais_lattice_type = str_sp.bravais_lattice_type;
  this->bravais_lattice_variation_type = str_sp.bravais_lattice_variation_type;
  this->bravais_lattice_system = str_sp.bravais_lattice_system;
  this->bravais_lattice_lattice_type = str_sp.bravais_lattice_lattice_type;
  this->bravais_lattice_lattice_variation_type = str_sp.bravais_lattice_lattice_variation_type;
  this->bravais_lattice_lattice_system = str_sp.bravais_lattice_lattice_system;
  this->volume_changed_original2new = str_sp.volume_changed_original2new; // DX20181024
  this->transform_coordinates_original2new = str_sp.transform_coordinates_original2new; // DX20181024
  this->transform_coordinates_new2original = str_sp.transform_coordinates_new2original; // DX20181024
  this->rotate_lattice_original2new = str_sp.rotate_lattice_original2new; // DX20181024
  this->rotate_lattice_new2original = str_sp.rotate_lattice_new2original; // DX20181024
  this->pearson_symbol = str_sp.pearson_symbol;
  this->crystal_family = str_sp.crystal_family;
  this->crystal_system = str_sp.crystal_system;
  this->point_group_crystal_class = str_sp.point_group_crystal_class;
  this->point_group_Shoenflies = str_sp.point_group_Shoenflies;
  this->point_group_Hermann_Mauguin = str_sp.point_group_Hermann_Mauguin;
  this->point_group_orbifold = str_sp.point_group_orbifold;
  this->point_group_type = str_sp.point_group_type;
  this->point_group_order = str_sp.point_group_order;
  this->point_group_structure = str_sp.point_group_structure;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [4] DONE" << endl;
  }

  // update sym_eps
  this->sym_eps = str_in.sym_eps = str_sc.sym_eps = str_sp.sym_eps; // DX
  this->sym_eps_calculated = str_in.sym_eps_calculated = str_sc.sym_eps_calculated = str_sp.sym_eps_calculated; // DX
  this->sym_eps_change_count = str_in.sym_eps_change_count = str_sc.sym_eps_change_count = str_sp.sym_eps_change_count;
  // DX20180222 - added sym_eps change count
  this->sym_eps_no_scan = str_in.sym_eps_no_scan = str_sc.sym_eps_no_scan = str_sp.sym_eps_no_scan;
  // DX20210430 - added no_scan
}

// ***************************************************************************
// Function GetReciprocalLatticeType //DX20210209
// ***************************************************************************
// Determine the reciprocal space symmetry information
// Includes self-consistency loop to ensure descriptions are commensurate
void xstructure::GetReciprocalLatticeType(double sym_eps) {
  xstructure str_sp;
  xstructure str_sc;
  GetReciprocalLatticeType(str_sp, str_sc, sym_eps);
}

void xstructure::GetReciprocalLatticeType(xstructure& str_sp, xstructure& str_sc, double sym_eps) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if (tolerance == AUROSTD_MAX_DOUBLE) {
    if ((*this).sym_eps_calculated) {
      tolerance = (*this).sym_eps;
    } else {
      tolerance = SYM::defaultTolerance((*this));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl;
  }

  // ---------------------------------------------------------------------------
  // RECIPROCAL - use klattice an one atom (at the origin)
  xstructure str_in;
  str_in.lattice = this->klattice;
  str_in.FixLattices();
  const _atom atom;
  str_in.AddAtom(atom, false); // CO20230319 - add by type

  // update tolerance info for all xstructures
  str_in.sym_eps = str_sp.sym_eps = str_sc.sym_eps = tolerance;
  str_in.sym_eps_calculated = str_sp.sym_eps_calculated = str_sc.sym_eps_calculated = (*this).sym_eps_calculated;
  str_in.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = (*this).sym_eps_change_count;
  str_in.sym_eps_no_scan = str_sp.sym_eps_no_scan = str_sc.sym_eps_no_scan = (*this).sym_eps_no_scan;
  // DX20210430 - added no_scan

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [1]" << endl;
  }
  // DX20170814 START - Use real pgroup to calculate pgroupk and then set pgroupk from str_sp to the pgroup and pgroup_xtal of str_reciprocal_in
  // DX20170814 The pgroup and pgroup_xtal are the same for the str_reciprocal structure because there is only one atom at the origin
  // DX20170814 (i.e. lattice and crystal symmetry are the same for the reciprocal space crystal)
  // DX20180426 - possible that lattice exhibits lower symmetry than crystal (i.e., from str_sp); would need to pass lattice symmetry from Standard_Lattice, but that information is not stored out of scope,
  // commenting out 5 lines below DX20170814 END
  LATTICE::Standard_Lattice_StructureDefault(str_in, str_sp, str_sc, false);
  // DX //DX20180226 - do not need to do full sym for recip

  this->reciprocal_lattice_type = str_sp.bravais_lattice_type;
  this->reciprocal_lattice_variation_type = str_sp.bravais_lattice_variation_type;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [2]" << endl;
  }

  // update sym_eps
  this->sym_eps = str_in.sym_eps = str_sc.sym_eps = str_sp.sym_eps; // DX
  this->sym_eps_calculated = str_in.sym_eps_calculated = str_sc.sym_eps_calculated = str_sp.sym_eps_calculated; // DX
  this->sym_eps_change_count = str_in.sym_eps_change_count = str_sc.sym_eps_change_count = str_sp.sym_eps_change_count;
  // DX20180222 - added sym_eps change count
  this->sym_eps_no_scan = str_in.sym_eps_no_scan = str_sc.sym_eps_no_scan = str_sp.sym_eps_no_scan;
  // DX20210430 - added no_scan
}

// ***************************************************************************
// Function GetSuperlatticeType() //DX20210302
// ***************************************************************************
// Determine the superlattice symmetry information
// Includes self-consistency loop to ensure descriptions are commensurate
void xstructure::GetSuperlatticeType(double sym_eps) {
  xstructure str_sp;
  xstructure str_sc;
  GetSuperlatticeType(str_sp, str_sc, sym_eps);
}

void xstructure::GetSuperlatticeType(xstructure& str_sp, xstructure& str_sc, double sym_eps) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  // ---------------------------------------------------------------------------
  // set symmetry tolerance based on the following sequence
  // 1) use input, 2) use sym_eps in xstructure, 3) calculate default
  double tolerance = sym_eps;
  if (tolerance == AUROSTD_MAX_DOUBLE) {
    if ((*this).sym_eps_calculated) {
      tolerance = (*this).sym_eps;
    } else {
      tolerance = SYM::defaultTolerance((*this));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [1] Set symmetry tolerance (starting sym_eps=" << tolerance << ")" << endl;
  }

  // ---------------------------------------------------------------------------
  // SUPERLATTICE - decorate with single atom type
  xstructure str_in = *this;
  str_in.ClearSymmetry(); // need to clear symmetry; otherwise, nothing is calculated
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [1]" << endl;
    cerr << str_in << endl;
  }
  // decorate with single atom type
  str_in.IdenticalAtoms(); // make superlattice
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [2]" << endl;
    cerr << str_in << endl;
  }
  // primitivize
  str_in.GetPrimitive(); // DX20210430 - remove obsolete eps=0.005
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [3]" << endl;
    cerr << str_in << endl;
  }
  // Minkowski
  str_in.Minkowski_calculated = false;
  str_in.MinkowskiBasisReduction();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [4]" << endl;
    cerr << str_in << endl;
  }

  // update tolerance info for all xstructures
  str_in.sym_eps = str_sp.sym_eps = str_sc.sym_eps = tolerance;
  str_in.sym_eps_calculated = str_sp.sym_eps_calculated = str_sc.sym_eps_calculated = (*this).sym_eps_calculated;
  str_in.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = (*this).sym_eps_change_count;
  str_in.sym_eps_no_scan = str_sp.sym_eps_no_scan = str_sc.sym_eps_no_scan = (*this).sym_eps_no_scan;
  // DX20210430 - added no_scan

  // main lattice function
  LATTICE::Standard_Lattice_StructureDefault(str_in, str_sp, str_sc, false);
  // DX //DX20180226 - do not need to do full sym for superlattice
  str_sp.ReScale(1.0); // DX20210211 - need to rescale to 1 since we aren't propagating the superlattice scaling factor
  this->bravais_superlattice_lattice = str_sp.lattice;
  this->bravais_superlattice_type = str_sp.bravais_lattice_type;
  this->bravais_superlattice_variation_type = str_sp.bravais_lattice_variation_type;
  this->bravais_superlattice_system = str_sp.bravais_lattice_system;
  this->pearson_symbol_superlattice = str_sp.pearson_symbol;

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [6] DONE" << endl;
  }

  // update sym_eps
  this->sym_eps = str_in.sym_eps = str_sc.sym_eps = str_sp.sym_eps; // DX
  this->sym_eps_calculated = str_in.sym_eps_calculated = str_sc.sym_eps_calculated = str_sp.sym_eps_calculated; // DX
  this->sym_eps_change_count = str_in.sym_eps_change_count = str_sc.sym_eps_change_count = str_sp.sym_eps_change_count;
  // DX20180222 - added sym_eps change count
  this->sym_eps_no_scan = str_in.sym_eps_no_scan = str_sc.sym_eps_no_scan = str_sp.sym_eps_no_scan;
  // DX20210430 - added no_scan
}

string GetLatticeType(xmatrix<double> lattice) {
  // DX double eps=0.00010,epsang=0.01;
  xstructure str_in;
  xstructure str_sp;
  xstructure str_sc;
  str_in.lattice = lattice;
  str_in.FixLattices();
  str_in.title = "NO_RECURSION";
  const _atom atom;
  str_in.AddAtom(atom, false); // CO20230319 - add by type
  // LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang); //SC OLD VERSION
  // DX int ss=0; //JX
  // DX LATTICE::Standard_Lattice_Structure(str_in,str_sp,str_sc,eps,epsang,ss,_EPS_); //JX
  LATTICE::Standard_Lattice_StructureDefault(str_in, str_sp, str_sc); // DX
  return str_sp.bravais_lattice_type;
  //  return str_sp.bravais_lattice_variation_type;
}

string GetLatticeType(xvector<double> data) {
  const xmatrix<double> lattice(GetClat(data));
  return GetLatticeType(lattice);
}

// ***************************************************************************
// Function Standard_Primitive
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure Standard_Primitive_UnitCellForm(const xstructure& a) {
  const xstructure str_in(a);
  xstructure str_sp;
  xstructure str_sc;
  if (str_in.Standard_Lattice_avoid == true) {
    return str_in; // Nothing to do
  }
  if (str_in.Standard_Lattice_primitive == true) {
    return str_in; // already Primitive
  }
  LATTICE::Standard_Lattice_StructureDefault(str_in, str_sp, str_sc);
  return str_sp;
}

xstructure GetStandardPrimitive(const xstructure& a) {
  // cerr << "GetStandardPrimitive(const xstructure& a)" << endl;
  return Standard_Primitive_UnitCellForm(a);
}

void xstructure::Standard_Primitive_UnitCellForm() {
  xstructure str_sp;
  xstructure str_sc;
  if (Standard_Lattice_avoid == true) {
    return; // Nothing to do
  }
  if (Standard_Lattice_primitive == true) {
    return; // already Primitive
  }
  // cerr << bravais_lattice_type << endl;
  LATTICE::Standard_Lattice_StructureDefault(*this, str_sp, str_sc);
  *this = str_sp;
}

void xstructure::GetStandardPrimitive() {
  Standard_Primitive_UnitCellForm();
}

xmatrix<double> GetStandardPrimitive(xmatrix<double> lattice) {
  xstructure str;
  str.lattice = lattice;
  str.scale = 1.0;
  str.FixLattices();
  str.title = "NO_RECURSION";
  const _atom atom;
  str.AddAtom(atom, false); // CO20230319 - add by type
  str = GetStandardPrimitive(str);
  return str.lattice;
}

xvector<double> GetStandardPrimitive(xvector<double> data) {
  const xmatrix<double> lattice(GetClat(data));
  return Getabc_angles(GetStandardPrimitive(lattice), DEGREES);
}

// ***************************************************************************
// Function Standard_Conventional
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure Standard_Conventional_UnitCellForm(const xstructure& a) {
  const xstructure str_in(a);
  xstructure str_sp;
  xstructure str_sc;
  if (str_in.Standard_Lattice_avoid == true) {
    return str_in; // Nothing to do
  }
  if (str_in.Standard_Lattice_conventional == true) {
    return str_in; // already Conventional
  }
  LATTICE::Standard_Lattice_StructureDefault(str_in, str_sp, str_sc);
  return str_sc;
}

xstructure GetStandardConventional(const xstructure& a) {
  return Standard_Conventional_UnitCellForm(a);
}

void xstructure::Standard_Conventional_UnitCellForm() {
  xstructure str_sp;
  xstructure str_sc;
  if (Standard_Lattice_avoid == true) {
    return; // Nothing to do
  }
  if (Standard_Lattice_conventional == true) {
    return; // already Conventional
  }
  LATTICE::Standard_Lattice_StructureDefault(*this, str_sp, str_sc);
  *this = str_sc;
}

void xstructure::GetStandardConventional() {
  Standard_Conventional_UnitCellForm();
}

xmatrix<double> GetStandardConventional(xmatrix<double> lattice) {
  xstructure str;
  str.lattice = lattice;
  str.scale = 1.0;
  str.FixLattices();
  str.title = "NO_RECURSION";
  const _atom atom;
  str.AddAtom(atom, false); // CO20230319 - add by type
  str = GetStandardConventional(str);
  return str.lattice;
}

xvector<double> GetStandardConventional(xvector<double> data) {
  const xmatrix<double> lattice(GetClat(data));
  return Getabc_angles(GetStandardConventional(lattice), DEGREES);
}

// ***************************************************************************
// Function SpeciesLabel
// ***************************************************************************
// Returns the name of the specie (the name of the 1st atom of the specie
// Stefano Curtarolo - nov 2008
string xstructure::SpeciesLabel(const uint& A) {
  if (A > num_each_type.size()) {
    return "xxx (outside boundaries)"; // outside boundaries
  }
  string name;
  vector<string> tokens;

  uint i = 0;
  uint A_start = 0; //,A_stop=0;
  for (i = 0; i < num_each_type.size(); i++) {
    if (i < A) {
      A_start += num_each_type[i];
    }
    // [UNNECESSARY] if(i==A) {A_start=A_start;}
    //    if(i==A) {A_stop=A_start+num_each_type[i];}
  }
  if (atoms.at(A_start).name_is_given == false) {
    return "xxx (name not given)"; // name not given
  }
  name = atoms.at(A_start).name;
  aurostd::StringSubstInPlace(name, "+", " ");
  aurostd::StringSubstInPlace(name, "-", " ");
  aurostd::string2tokens(name, tokens, " ");
  name = tokens[0];
  return name;
}

string SpeciesLabel(const xstructure& str, const uint& A) {
  xstructure sstr(str);
  return sstr.SpeciesLabel(A);
}

// ***************************************************************************
// Function SpeciesSwap
// ***************************************************************************
// Permute Species A with B (safe for species C) Stefano Curtarolo Nov 2008
void xstructure::SpeciesSwap(const uint& specieA, const uint& specieB) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  // some useful checks
  if (num_each_type.empty()) {
    return; // empty structures have nothing to swap
  }
  if (num_each_type.size() == 1) {
    return; // pure structures have nothing to swap
  }
  // check done
  deque<_atom> atoms_buf;
  uint i = 0;
  uint j = 0;
  uint k = 0;
  uint specieAA = 0;
  uint specieAA_start = 0;
  uint specieAA_stop = 0;
  uint specieBB = 0;
  uint specieBB_start = 0;
  uint specieBB_stop = 0;
  // PUT ALPHABETIC
  if (specieA < specieB) {
    specieAA = specieA;
    specieBB = specieB;
  }
  if (specieA > specieB) {
    specieAA = specieB;
    specieBB = specieA;
  }
  // CHECK
  if (specieA == specieB) {
    return; // NO SWAP
  }
  if (specieA > num_each_type.size()) {
    return; // NO SWAP nothing to swap
  }
  if (specieB > num_each_type.size()) {
    return; // NO SWAP nothing to swap
  }
  // FIND BOUNDARIES
  for (i = 0; i < num_each_type.size(); i++) {
    if (i < specieAA) {
      specieAA_start += num_each_type[i];
    }
    if (i == specieAA) {
      specieAA_stop = specieAA_start + num_each_type[i];
    }
    if (i < specieBB) {
      specieBB_start += num_each_type[i];
    }
    if (i == specieBB) {
      specieBB_stop = specieBB_start + num_each_type[i];
    }
  }
  if (LDEBUG) {
    cerr << "DEBUG  atoms.size()=" << atoms.size() << endl;
  }
  if (LDEBUG) {
    cerr << "DEBUG  specieAA=" << specieAA << endl;
  }
  if (LDEBUG) {
    cerr << "DEBUG  specieBB=" << specieBB << endl;
  }
  if (LDEBUG) {
    cerr << "DEBUG  specieAA_start=" << specieAA_start << endl;
  }
  if (LDEBUG) {
    cerr << "DEBUG  specieAA_stop=" << specieAA_stop << endl;
  }
  if (LDEBUG) {
    cerr << "DEBUG  specieBB_start=" << specieBB_start << endl;
  }
  if (LDEBUG) {
    cerr << "DEBUG  specieBB_stop=" << specieBB_stop << endl;
  }
  // ATOMS  -----------------------------
  atoms_buf.clear();
  for (i = 0; i < atoms.size(); i++) {
    atoms_buf.push_back(atoms[i]); // create buffer
  }
  atoms.clear(); // RECONSTRUCTING
  for (i = 0; i < specieAA_start; i++) {
    atoms.push_back(atoms_buf.at(i)); // before specieA preserve
  }
  for (i = specieBB_start; i < specieBB_stop; i++) {
    atoms.push_back(atoms_buf.at(i)); // during specieA swap with specieB
  }
  for (i = specieAA_stop; i < specieBB_start; i++) {
    atoms.push_back(atoms_buf.at(i));
  }
  // between specieA and specieB preserve
  for (i = specieAA_start; i < specieAA_stop; i++) {
    atoms.push_back(atoms_buf.at(i)); // during specieB swap with specieA
  }
  for (i = specieBB_stop; i < atoms_buf.size(); i++) {
    atoms.push_back(atoms_buf.at(i)); // after specieB preserve
  }
  // swap numbers
  uint iaus;
  iaus = num_each_type.at(specieAA);
  num_each_type.at(specieAA) = num_each_type.at(specieBB);
  num_each_type.at(specieBB) = iaus;
  // now fix the types
  k = 0;
  for (i = 0; i < num_each_type.size(); i++) {
    for (j = 0; j < (uint) num_each_type[i]; j++) {
      atoms.at(k++).type = i; // done
    }
  }
  // swap species
  string saus;
  double daus;
  if (specieAA < species.size() && specieBB < species.size()) {
    saus = species.at(specieAA);
    species.at(specieAA) = species.at(specieBB);
    species.at(specieBB) = saus;
  }
  if (specieAA < species_pp.size() && specieBB < species_pp.size()) {
    saus = species_pp.at(specieAA);
    species_pp.at(specieAA) = species_pp.at(specieBB);
    species_pp.at(specieBB) = saus;
  }
  if (specieAA < comp_each_type.size() && specieBB < comp_each_type.size()) {
    daus = comp_each_type.at(specieAA);
    comp_each_type.at(specieAA) = comp_each_type.at(specieBB);
    comp_each_type.at(specieBB) = daus;
  } // CO20180705
  if (specieAA < stoich_each_type.size() && specieBB < stoich_each_type.size()) {
    daus = stoich_each_type.at(specieAA);
    stoich_each_type.at(specieAA) = stoich_each_type.at(specieBB);
    stoich_each_type.at(specieBB) = daus;
  } // CO20180705
  if (specieAA < species_volume.size() && specieBB < species_volume.size()) {
    daus = species_volume.at(specieAA);
    species_volume.at(specieAA) = species_volume.at(specieBB);
    species_volume.at(specieBB) = daus;
  }
  if (specieAA < species_mass.size() && specieBB < species_mass.size()) {
    daus = species_mass.at(specieAA);
    species_mass.at(specieAA) = species_mass.at(specieBB);
    species_mass.at(specieBB) = daus;
  }
}

// ***************************************************************************
// Function SpeciesGetAlphabetic
// ***************************************************************************
// Tell if two species are alphabetic!  Stefano Curtarolo Nov 2008
bool xstructure::SpeciesGetAlphabetic() {
  stringstream message;
  if (num_each_type.size() != species.size()) {
    message << "num_each_type.size()!=species.size()   (" << num_each_type.size() << "," << species.size() << ")" << endl;
    message << "num_each_type.size()=" << num_each_type.size() << ": ";
    for (size_t i = 0; i < num_each_type.size(); i++) {
      message << num_each_type[i] << " ";
    }
    message << endl;
    message << "species.size()=" << species.size() << ": ";
    for (size_t i = 0; i < species.size(); i++) {
      message << species[i] << " ";
    }
    message << endl;
    message << "species_pp.size()=" << species_pp.size() << ": ";
    for (size_t i = 0; i < species_pp.size(); i++) {
      message << species_pp[i] << " ";
    }
    message << endl;
    message << "species_volume.size()=" << species_volume.size() << ": ";
    for (size_t i = 0; i < species_volume.size(); i++) {
      message << species_volume[i] << " ";
    }
    message << endl;
    message << "species_mass.size()=" << species_mass.size() << ": ";
    for (size_t i = 0; i < species_mass.size(); i++) {
      message << species_mass[i] << " ";
    }
    message << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_RANGE_);
  }
  // some useful checks
  if (species.empty()) {
    return true; // empty structures are always alphabetic
  }
  if (species.size() == 1) {
    return true; // pure structures are always alphabetic
  }
  // check done
  vector<string> sspecies;
  for (size_t isp = 0; isp < species.size(); isp++) {
    sspecies.push_back(aurostd::RemoveNumbers(KBIN::VASP_PseudoPotential_CleanName(species[isp])));
  }

  for (size_t isp = 0; isp < sspecies.size() - 1; isp++) {
    if (sspecies[isp] > sspecies.at(isp + 1)) {
      return false;
    }
  }
  // otherwise return true;
  return true;
}

// ***************************************************************************
// Function SpeciesPutAlphabetic
// ***************************************************************************
// Tell if two species are alphabetic!  Stefano Curtarolo Nov 2008
bool xstructure::SpeciesPutAlphabetic() {
  stringstream message;
  if (num_each_type.size() != species.size()) {
    message << "num_each_type.size()!=species.size()   (" << num_each_type.size() << "," << species.size() << ")";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_RANGE_);
  }
  // some useful checks
  if (species.empty()) {
    return true; // empty structures are always alphabetic
  }
  if (species.size() == 1) {
    return true; // pure structures are always alphabetic
  }
  // check done
  vector<string> sspecies;

  while (SpeciesGetAlphabetic() == false) {
    sspecies.clear();
    for (size_t isp = 0; isp < species.size(); isp++) {
      sspecies.push_back(aurostd::RemoveNumbers(KBIN::VASP_PseudoPotential_CleanName(species[isp])));
    }
    for (size_t isp = 0; isp < sspecies.size() - 1; isp++) {
      if (sspecies[isp] > sspecies.at(isp + 1)) {
        SpeciesSwap(isp, isp + 1);
      }
    }
  }
  MakeBasis(); // repetita iuvant
  // otherwise return true;
  return true;
}

// ***************************************************************************
// Function SpeciesString
// ***************************************************************************
// Returns a string with the species.  Stefano Curtarolo Nov 2008
string xstructure::SpeciesString() {
  stringstream strstream;
  strstream.clear();
  strstream.str(std::string());
  for (size_t i = 0; i < num_each_type.size(); i++) {
    strstream << SpeciesLabel(i);
    if (i != num_each_type.size() - 1) {
      strstream << " ";
    }
  }
  return strstream.str();
}

// ***************************************************************************
// Function SetSpecies
// ***************************************************************************
// Set the species  Stefano Curtarolo Nov 2014
// CO20221112 - do not set species explicitly, always use add/remove atoms, it updates EVERYTHING all at once and avoids seg faults
uint xstructure::SetSpecies(const deque<string>& species_in, bool sort_species) { // CO20230319 - sort species
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_in=" << aurostd::joinWDelimiter(species_in, ",") << endl;
    cerr << __AFLOW_FUNC__ << " xstr=" << endl << (*this) << endl;
    cerr << __AFLOW_FUNC__ << " xstr.species=" << aurostd::joinWDelimiter(species, ",") << endl;
    cerr << __AFLOW_FUNC__ << " xstr.num_each_type=" << aurostd::joinWDelimiter(num_each_type, ",") << endl;
  }
  if (species_in.size() != species.size()) {
    aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "species_in.size()!=species.size()", _VALUE_RANGE_);
  } // CO20190317
  if (species_in.size() != num_each_type.size()) {
    aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "species_in.size()!=num_each_type.size()", _VALUE_RANGE_);
  } // CO20190317
  uint iatom = 0;
  deque<_atom> atoms_new;
  for (size_t itype = 0; itype < num_each_type.size(); itype++) {
    for (int j = 0; j < num_each_type[itype]; j++) {
      atoms_new.push_back(atoms[iatom++]);
      atoms_new.back().name = species_in[itype];
      atoms_new.back().CleanName();
      //[DX20170921 - Need to keep spin info]atoms_new.back().CleanSpin();
      atoms_new.back().name_is_given = (!atoms_new.back().name.empty());
    }
  }
  ReplaceAtoms(atoms_new, true, sort_species); // CO20230319 - sort species
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " xstr_updated=" << endl << (*this) << endl;
    cerr << __AFLOW_FUNC__ << " xstr_updated.species=" << aurostd::joinWDelimiter(species, ",") << endl;
    cerr << __AFLOW_FUNC__ << " xstr_updated.num_each_type=" << aurostd::joinWDelimiter(num_each_type, ",") << endl;
  }
  return species_in.size();
}

// ***************************************************************************
// Function UpdateSpecies() //DX20210202 [from AddAtom, consolidate to function]
// ***************************************************************************
void xstructure::UpdateSpecies(const _atom& atom, bool add_species) { // CO20230319 - fixing species vs. types issue

  // Update the species info based on the atom input
  // If the species is already in xstructure, update the number of types
  // and composition of each type, otherwise, add the new species info
  // This code was copied from AddAtom

  const bool LDEBUG = (false || XHOST.DEBUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    cerr << __AFLOW_FUNC__ << " add_species=" << add_species << endl;
  }

  bool FOUND_SPECIES = false;
  uint species_position = 0;

  if (add_species) { // CO20230319 - if we want to decide where to add the atom based on atom.name
    if (atom.name.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "atom.name.empty()", _INPUT_MISSING_);
    } // DX20210324 - check if name is empty
    for (size_t isp = 0; isp < species.size() && FOUND_SPECIES == false; isp++) {
      if (atom.name == species[isp]) {
        FOUND_SPECIES = true;
        species_position = isp;
      }
    }
    if (!FOUND_SPECIES) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " new_species=" << atom.name << endl;
      }
      species_position = species.size();
    }
  } else { // DX20210324 - if name is empty, then use types
    for (size_t isp = 0; isp < species.size() && FOUND_SPECIES == false; isp++) {
      if (atom.type == (int) isp) {
        FOUND_SPECIES = true;
        species_position = isp;
      }
    }
    if (!FOUND_SPECIES) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " new_type=" << atom.type << endl;
      }
      species_position = atom.type;
    }
  }

  if (FOUND_SPECIES == false) {
    while (species.size() <= species_position) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " increasing species.size()=" << species.size() << endl;
      }
      num_each_type.push_back(0);
      comp_each_type.push_back(0.0);
      species.emplace_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp.emplace_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp_type.emplace_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp_version.emplace_back(""); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp_ZVAL.push_back(0.0); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_pp_vLDAU.emplace_back(); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_volume.push_back(GetAtomVolume("")); // cerr << "UpdateSpecies=" << atom.name << endl;
      species_mass.push_back(GetAtomMass("")); // cerr << "UpdateSpecies=" << atom.name << endl;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " species.size()=" << species.size() << endl;
      }
    }
  }
  //
  if (num_each_type[species_position] == 0) {
    species[species_position] = atom.name; // cerr << "UpdateSpecies=" << atom.name << endl;
    species_pp[species_position] = atom.name; // cerr << "UpdateSpecies=" << atom.name << endl;
    species_volume[species_position] = GetAtomVolume(atom.name); // cerr << "UpdateSpecies=" << atom.name << endl;
    species_mass[species_position] = GetAtomMass(atom.name); // cerr << "UpdateSpecies=" << atom.name << endl;
  }
  num_each_type[species_position]++;
  comp_each_type[species_position] += atom.partial_occupation_value;
}

// ***************************************************************************
// NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGLI NIGGL
// ***************************************************************************
// GetNiggliCell
// ***************************************************************************
// Calculates the reduced Niggli cell.  It is based on
// a program of Eric Wu's.  Here is his + my documentation.
// FUNCTION:
// Subroutine calculates the reduced cell (Niggli cell) for a given
// primitive cell.  It also gives the transformation matrix to go from
// the primitive cell to the Niggli cell.
//
// REFERENCES:
// Y. Lepage, J. Appl. Cryst 15 255-229 (1982) and
// Y. Lepage, J. Appl. Cryst 20 264-269 (1987)
// International tables for crystallography Ch 5,9
// I Krivy, B. Gruber. Acta Cryst 1976 A32 297 (1975)
// W. Clegg., Acta Cryst A37 913-915 (1981)
// A.D. Mighell J. Appl Cryst 9 491-498 (1976)
//
// METHODOLOGY:
// The Niggli cell is UNIQUE and reduced. A Niggli cell is a Buerger cell
//(shortest lattice vectors) but  not necessarily the other way around,
// since the Buerger cell is not unique.  This info is useful because with
// the Niggli cell one knows instantly the type of bravais lattice
// and how to go to the conventional cell if you apply the correct
// algorithm.  There are 2 ways to do this.  The first is the straight
// way -  see International Tables Ch. 9 .  The other involves counting
// the # of 2fold axes and going from there, see the Lepage and Clegg
// references.  This second method is more useful if you want to find
// more than just bravais lattice (symmetry elements or space group
// for example).  The subroutine assumes you input a primitive cell
// and it spits out the Niggli cell and the transformation matrix
// to go from the primitive cell to Niggli cell.
// NB! There is a typo in step 3 of Krivy and Gruber that I corrected.
//[ (ksi*eta*ksi) should be (ksi*eta*zeta) in first part of Eq. 3. ]
//
// INPUT:
// in_lat: a primitive set of lattice vectors,
// in_lat(1,1)=a_x
// in_lat(1,2)=a_y
// in_lat(2,1)=b_x  ...
//
// OUTPUT:
// niggli_lat - reduced (Niggli) cell, same format
// P - Matrix to go from Primitive to Niggli cell
// Q - P^-1;
//
// For more info of the definitions here see Int. Tables Cryst. ch.5
// Let nlat be the niggli lattice, where the first column is a,
// second b, third c.  Let lat be the original cell equivalent.
// The P matrix to transform cells is generally defined by
// nlat = lat * P
// We will follow this convention for defining P and Q.
// Note that is the rest of the code, we define lattices with
// rows equal to a, b, and c.  Therefore, under our conventional
// definitions of nlat, call it nlat_c=nlat^T, and lat, call
// it lat_c=lat^T, we have
// nlat_c = P^T * lat_c
// where ^T means the transpose.

bool GetNiggliCell(const xmatrix<double>& in_lat, xmatrix<double>& niggli_lat, xmatrix<double>& P, xmatrix<double>& Q) {
  // return false if failed
  // rerurn true if ok

  //  double RENORM=1E+6;  //DM not used
  const int MAXITER = 100000;
  const double TOL = _ZERO_TOL_; // DX20180212 - from 1e-15 to 1e-10

  // Initialize matrices for tranformations (3x3).
  const xmatrix<double> m1(3, 3);
  m1(1, 2) = 1.0;
  m1(2, 1) = 1.0;
  m1(3, 3) = -1.0;
  const xmatrix<double> m2(3, 3);
  m2(1, 1) = -1.0;
  m2(2, 3) = 1.0;
  m2(3, 2) = 1.0;
  const xmatrix<double> m3(3, 3);
  // Els 1,1 2,2 3,3 are changed later for m3
  const xmatrix<double> m4(3, 3);
  // Els 1,1 2,2 3,3 are changed later for m4
  const xmatrix<double> m5(3, 3);
  m5(1, 1) = 1.0;
  m5(2, 2) = 1.0;
  m5(3, 3) = 1.0;
  // El 3,2 is changed later for m5
  const xmatrix<double> m6(3, 3);
  m6(1, 1) = 1.0;
  m6(2, 2) = 1.0;
  m6(3, 3) = 1.0;
  // El 3,1 is changed later for m6
  const xmatrix<double> m7(3, 3);
  m7(1, 1) = 1.0;
  m7(2, 2) = 1.0;
  m7(3, 3) = 1.0;
  // El 2,1 is changed later for m7
  const xmatrix<double> m8(3, 3);
  m8(1, 1) = 1.0;
  m8(2, 2) = 1.0;
  m8(1, 3) = 1.0;
  m8(2, 3) = 1.0;
  m8(3, 3) = 1.0;

  // Initialize a, b, c, ksi, eta, zeta
  xvector<double> indat(6);
  indat = Getabc_angles(in_lat, RADIANS);
  double a = indat(1) * indat(1);
  double b = indat(2) * indat(2);
  double c = indat(3) * indat(3);
  double ksi = 2.0 * indat(2) * indat(3) * cos(indat(4));
  double eta = 2.0 * indat(1) * indat(3) * cos(indat(5));
  double zeta = 2.0 * indat(1) * indat(2) * cos(indat(6));
  // Round to RENORM decimal place to eliminate numerical errors
  //  double a=Nint(indat(1)*indat(1)*RENORM);
  //  double b=Nint(indat(2)*indat(2)*RENORM);
  //  double c=Nint(indat(3)*indat(3)*RENORM);
  //  double ksi=Nint(2.0*indat(2)*indat(3)*cos(indat(4))*RENORM);
  //  double eta=Nint(2.0*indat(1)*indat(3)*cos(indat(5))*RENORM);
  //  double zeta=Nint(2.0*indat(1)*indat(2)*cos(indat(6))*RENORM);

  // Dummy variables
  double temp;
  double temp1;
  double temp2;
  double temp3;

  // Initialize tranformation matrix
  //  xmatrix<double> P(3,3);
  P(1, 1) = 1;
  P(2, 2) = 1;
  P(3, 3) = 1;

  int cnt = 0;

  // Start loop

LoopHead:

  cnt++;
  //  cerr << "DEBUG: Niggli() cnt=" << cnt << endl;
  if (cnt > MAXITER) {
    //     stringstream oss;
    //     // -------
    //     oss << "EEEEE  aflow " << VERSION << endl;
    //     oss << "EEEEE  ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    //     oss << "EEEEE  ERROR: Number of interations greater the MAXITER = " << MAXITER << endl;
    //     oss << "EEEEE  ERROR: This seems like too many - there is probably some problem." << endl;
    //     oss << "EEEEE  ERROR:." << endl;
    //     //  oss << endl;
    //     // -------
    //     cerr << oss.str();
    //     //  cout << oss.str();
    //     // -------
    return false;
  }

  // Step 1
  if (((a - b) > TOL) || ((std::abs(a - b) < TOL) && ((std::abs(ksi) - std::abs(eta)) > TOL))) // DX20180209 - more robust; precision
  { // CO20200106 - patching for auto-indenting
    temp = a;
    a = b;
    b = temp;
    temp = -ksi;
    ksi = -eta;
    eta = temp;
    P = P * m1;
  }

  // Step 2
  if (((b - c) > TOL) || ((std::abs(b - c) < TOL) && ((std::abs(eta) - std::abs(zeta)) > TOL))) // DX20180209 - more robust; precision
  { // CO20200106 - patching for auto-indenting
    temp = b;
    b = c;
    c = temp;
    temp = -eta;
    eta = -zeta;
    zeta = temp;
    P = P * m2;
    goto LoopHead;
  }

  // Step 3
  if ((ksi * eta * zeta) > TOL) // DX20180209 - more robust; precision
  { // CO20200106 - patching for auto-indenting
    m3(1, 1) = SignNoZero(ksi);
    m3(2, 2) = SignNoZero(eta);
    m3(3, 3) = SignNoZero(zeta);
    P = P * m3;
    ksi = std::abs(ksi);
    eta = std::abs(eta);
    zeta = std::abs(zeta);
  }

  // Step 4
  // EG's original if((abs(ksi)<TOL)||(abs(eta)<TOL)||(abs(zeta)<TOL)||(ksi*eta*zeta<=0))
  if ((std::abs(ksi) < TOL) || (std::abs(eta) < TOL) || (std::abs(zeta) < TOL) || (ksi * eta * zeta <= -TOL))
  // DX20180212 - if any are zero, then, ksi*eta*zeta is zero
  // if(ksi*eta*zeta<=TOL) //DX20180209 - more robust; precision
  { // CO20200106 - patching for auto-indenting
    m4(1, 1) = -(SignNoZero(ksi));
    m4(2, 2) = -(SignNoZero(eta));
    m4(3, 3) = -(SignNoZero(zeta));
    if (std::abs(ksi) < TOL) {
      m4(1, 1) = m4(2, 2) * m4(3, 3);
    }
    if (std::abs(eta) < TOL) {
      m4(2, 2) = m4(1, 1) * m4(3, 3);
    }
    if (std::abs(zeta) < TOL) {
      m4(3, 3) = m4(1, 1) * m4(2, 2);
    }
    P = P * m4;
    ksi = -std::abs(ksi);
    eta = -std::abs(eta);
    zeta = -std::abs(zeta);
  }

  // Step 5
  if (((std::abs(ksi) - b) > TOL) || ((std::abs(ksi - b) < TOL) && (2.0 * eta - zeta) < -TOL) || // DX20180209 - more robust; precision
      ((std::abs(ksi + b) < TOL) && (zeta < -TOL))) // DX20180209 - more robust; precision
  { // CO20200106 - patching for auto-indenting
    m5(2, 3) = -(SignNoZero(ksi));
    P = P * m5;
    temp1 = b + c - ksi * SignNoZero(ksi);
    temp2 = eta - zeta * SignNoZero(ksi);
    temp3 = ksi - 2.0 * b * SignNoZero(ksi);
    c = temp1;
    eta = temp2;
    ksi = temp3;
    goto LoopHead;
  }

  // Step 6
  if (((std::abs(eta) - a) > TOL) || ((std::abs(eta - a) < TOL) && (2.0 * ksi - zeta) < -TOL) || // DX20180209 - more robust; precision
      ((std::abs(eta + a) < TOL) && (zeta < -TOL))) // DX20180209 - more robust; precision
  { // CO20200106 - patching for auto-indenting
    m6(1, 3) = -SignNoZero(eta);
    P = P * m6;
    temp1 = a + c - eta * SignNoZero(eta);
    temp2 = ksi - zeta * SignNoZero(eta);
    temp3 = eta - 2.0 * a * SignNoZero(eta);
    c = temp1;
    ksi = temp2;
    eta = temp3;
    goto LoopHead; // DX20180212 - this was missing
  }

  // Step 7
  if (((std::abs(zeta) - a) > TOL) || ((std::abs(zeta - a) < TOL) && (2.0 * ksi - eta) < -TOL) || // DX20180209 - more robust; precision
      ((std::abs(zeta + a) < TOL) && (eta < -TOL))) // DX20180209 - more robust; precision
  { // CO20200106 - patching for auto-indenting
    m7(1, 2) = -SignNoZero(zeta);
    P = P * m7;
    temp1 = a + b - zeta * SignNoZero(zeta);
    temp2 = ksi - eta * SignNoZero(zeta);
    temp3 = zeta - 2.0 * a * SignNoZero(zeta);
    b = temp1;
    ksi = temp2;
    zeta = temp3;
    goto LoopHead;
  }

  // Step 8
  if ((ksi + eta + zeta + a + b < -TOL) || // DX20180209 - more robust; precision
      (std::abs(ksi + eta + zeta + a + b) < TOL && (2.0 * (a + eta) + zeta > TOL))) // DX20180209 - more robust; precision
  { // CO20200106 - patching for auto-indenting
    P = P * m8;
    temp1 = a + b + c + ksi + eta + zeta;
    temp2 = 2.0 * b + ksi + zeta;
    temp3 = 2.0 * a + eta + zeta;
    c = temp1;
    ksi = temp2;
    eta = temp3;
    goto LoopHead;
  }

  // Renormalize back to regular cell (divide by RENORM)
  const xvector<double> outdat(6);
  //      outdat(1)=sqrt(a/RENORM);
  //      outdat(2)=sqrt(b/RENORM);
  //      outdat(3)=sqrt(c/RENORM);
  //      outdat(4)=acos(ksi/RENORM/2.0/outdat(2)/outdat(3));
  //      outdat(5)=acos(eta/RENORM/2.0/outdat(1)/outdat(3));
  //      outdat(6)=acos(zeta/RENORM/2.0/outdat(1)/outdat(2));
  outdat(1) = sqrt(a);
  outdat(2) = sqrt(b);
  outdat(3) = sqrt(c);
  outdat(4) = acos(ksi / 2.0 / outdat(2) / outdat(3));
  outdat(5) = acos(eta / 2.0 / outdat(1) / outdat(3));
  outdat(6) = acos(zeta / 2.0 / outdat(1) / outdat(2));

  // Get Niggli cell
  niggli_lat = trasp(P) * in_lat;

  // Get Q
  //  Q = P;
  // xmatrix<double> tmat(3,3);
  //  GaussJordan(Q,tmat); // Returns Q=P^-1.
  Q = inverse(P);

  // Checks
  // Make sure that a,b,c,alpha,beta,gamma are the same from
  // direct calculation and from using P to get niggli_lat.
  xvector<double> poutdat(6);
  poutdat = Getabc_angles(niggli_lat, RADIANS);
  int flag = 0;
  for (int i = 1; i <= 6; i++) {
    if (std::abs(poutdat(i) - outdat(i)) > 2 * TOL) {
      flag = 1;
    }
  }
  if (flag) {
    stringstream cout;
    cout << "ERROR: CellReduceFuncs/GetNiggliCell" << endl;
    cout << "ERROR: Lattice parameters/angles as calculated" << endl;
    cout << "ERROR: with Niggli algorithm and P do not match." << endl;
    cout << "ERROR: a,b,c,alpha,beta,gamma from direct algorithm: ";
    cout << outdat << endl;
    cout << "ERROR: a,b,c,alpha,beta,gamma from P matrix: ";
    cout << poutdat << endl;
    cout << "ERROR: Returning." << endl;
    // output
    cerr << cout.str();
    //  cout << oss.str();
    // done
    return false;
  }
  return true; // perfect
}

// DX20180213 - Niggli algorithm, fixed tolerances and added goto loop in step 6 - END

xstructure GetNiggliStr(const xstructure& in_str) {
  xmatrix<double> niggli_lat(3, 3);
  xmatrix<double> P(3, 3);
  xmatrix<double> Q(3, 3);
  xstructure sstr = in_str;
  const double scale = sstr.scale;
  sstr = ReScale(sstr, 1.0);
  bool is_Niggli_ok = true;
  is_Niggli_ok = GetNiggliCell(sstr.lattice, niggli_lat, P, Q);
  if (is_Niggli_ok == false) {
    sstr = in_str;
    sstr.Niggli_has_failed = true;
    return sstr; // got a problem
  } else {
    // Create new str same as before but with Niggli cell params.
    // Transform cell parameters by xyz_niggli = Q * xyz_original
    // (see ITC, ch.5).
    xstructure niggli_str = sstr;
    niggli_str.Niggli_has_failed = false;
    niggli_str.lattice = niggli_lat;
    niggli_str.FixLattices();
    niggli_str.atoms.clear();

    _atom atom;
    for (int ia = 0; ia < (int) sstr.atoms.size(); ia++) {
      atom = sstr.atoms[ia];
      atom.fpos = Q * sstr.atoms[ia].fpos;
      atom.cpos = F2C(niggli_str.lattice, atom.fpos);
      niggli_str.atoms.push_back(atom);
    }
    // Set lattice params, atom positions, put in the 0th cell, and reset scale.
    niggli_str = BringInCell(niggli_str);
    for (int ia = 0; ia < (int) sstr.atoms.size(); ia++) { // NEED TO CLEAR lattice positions of atoms
      clear(sstr.atoms[ia].ijk); // because the niggly reproduces the same structure
    }
    // DONE and fix it back
    niggli_str = ReScale(niggli_str, scale);
    // DONE
    niggli_str.Niggli_calculated = true;
    return niggli_str;
  }
}

xmatrix<double> GetNiggliStr(const xmatrix<double>& lattice) {
  xmatrix<double> niggli_lat(3, 3);
  xmatrix<double> P(3, 3);
  xmatrix<double> Q(3, 3);
  bool is_Niggli_ok = true;
  is_Niggli_ok = GetNiggliCell(lattice, niggli_lat, P, Q);
  if (is_Niggli_ok == false) {
    return lattice;
  }
  return niggli_lat;
}

// ***************************************************************************
// Function NiggliUnitCellForm
// ***************************************************************************
// Converts the unit cell to the standardized Niggli form.
// It is based on the GetNiggliStr() function written by Eric Wu and
// Dane Morgan, adapted by Stefano Curtarolo and present in aflow_pflow_funcs.cpp
xstructure NiggliUnitCellForm(const xstructure& a) {
  xstructure b(a);
  b.NiggliUnitCellForm();
  return b;
}

void xstructure::NiggliUnitCellForm() {
  if (Niggli_avoid == true) {
    return;
  }
  if (Niggli_calculated == false) {
    *this = GetNiggliStr(*this);
    Niggli_calculated = true;
    FixLattices();
    Niggli_has_failed = false;
  } else {
    if (!XHOST.QUIET) {
      cout << XPID << "00000  MESSAGE NIGGLI Form already Calculated: skipping" << endl;
    }
  }
}

xmatrix<double> NiggliUnitCellForm(const xmatrix<double>& lattice) {
  return GetNiggliStr(lattice);
}

// ***************************************************************************
// Function GetNiggliStructures() //DX20201006
// ***************************************************************************
void GetNiggliStructures(vector<xstructure>& structures, uint start_index, uint end_index) {
  // Converts a set of xstructures to their Niggli representation
  // Optional indices can be included; useful for pre-distributed
  // threading schemes
  // Default: run over entire range

  // if end index is greater than structures.size(), then compute Niggli cell for all structures
  if (end_index > structures.size()) {
    end_index = structures.size();
  }

  for (uint i = start_index; i < end_index; i++) {
    structures[i].NiggliUnitCellForm();
  }
}

// **************************************************************************
// Function MinkowskiBasisReduction
// **************************************************************************
// This routine takes a set of basis vectors (that form a lattice)
// and reduces them so that they form the shortest possible basis.
// The reduction is performed so that each vector "a_i" is a close
// as possible to the origin while remaining in the affine plane which
// is defined by "a_j", "a_k" but shifted by "a_i", for any choice
// of even permutations of i,j,k in 1,2,3.
// See Lecture notes in computer science, ISSN 0302-974, ANTS - VI :
// algorithmic number theory, 2004, vol. 3076, pp. 338-357
// ISBN 3-540-22156-5
// Written by Gus Hart in F90, recoded by SC in C++ (Sep/08).
xstructure MinkowskiBasisReduction(const xstructure& a) {
  xstructure b(a);
  b.MinkowskiBasisReduction();
  return b;
}

void xstructure::MinkowskiBasisReduction() {
  if (LatticeReduction_avoid == true) {
    return;
  }
  if (Minkowski_avoid == true) {
    return;
  }
  if (Minkowski_calculated == false) {
    xmatrix<double> basis(3, 3);
    const double old_scale = scale;
    ReScale(1.0);
    basis = trasp(lattice);
    basis = aurostd::reduce_to_shortest_basis(basis);
    lattice = trasp(basis);
    FixLattices(); // get f2c c2f and klattice
    for (size_t i = 0; i < atoms.size(); i++) {
      atoms[i].fpos = C2F(lattice, atoms[i].cpos);
    }
    // atoms[i]=C2F(lattice,atoms[i]);  // it works, same just curiosity
    // atoms[i]=C2F(atoms[i]);          // it works, same just curiosity
    ReScale(old_scale);
    FixLattices();
    Minkowski_calculated = true;
    Minkowski_has_failed = false;
  } else {
    //    if(!QUIET) cout << XPID << "00000  MESSAGE MINKOWSKI Basis Reduction already Calculated: skipping" << endl;
  }
}

xmatrix<double> MinkowskiBasisReduction(const xmatrix<double>& lattice) {
  xmatrix<double> basis(3, 3);
  basis = lattice;
  basis = trasp(basis);
  basis = aurostd::reduce_to_shortest_basis(basis);
  basis = trasp(basis);
  return basis;
}

// ***************************************************************************
// Function GetMinkowskiStructures() //DX20201006
// ***************************************************************************
void GetMinkowskiStructures(vector<xstructure>& structures, uint start_index, uint end_index) {
  // Converts a set of xstructures to their Minkowski representation
  // Optional indices can be included; useful for pre-distributed
  // threading schemes
  // Default: run over entire range

  // if end index is greater than structures.size(), then compute Minkowski cell for all structures
  if (end_index > structures.size()) {
    end_index = structures.size();
  }

  for (uint i = start_index; i < end_index; i++) {
    structures[i].MinkowskiBasisReduction();
  }
}

// ***************************************************************************
// Function LatticeReduction
// ***************************************************************************
// Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
xstructure LatticeReduction(const xstructure& a) {
  xstructure b(a);
  b.LatticeReduction();
  return b;
}

void xstructure::LatticeReduction() {
  if (LatticeReduction_avoid == true) {
    return;
  }
  if (LatticeReduction_calculated == false) {
    MinkowskiBasisReduction(); // Minkowski first
    NiggliUnitCellForm(); // Niggli Second
    FixLattices(); // fix f2c c2f and K-space
    BringInCell(); // bring atoms in new basis
    LatticeReduction_calculated = true;
    LatticeReduction_has_failed = false;
  } else {
    //   if(!QUIET) cout << XPID << "00000  MESSAGE LATTICE Basis Reduction already Calculated: skipping" << endl;
  }
}

xmatrix<double> LatticeReduction(const xmatrix<double>& lattice) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  xmatrix<double> basis(3, 3);
  basis = lattice;
  basis = MinkowskiBasisReduction(basis); // Minkowski first
  basis = NiggliUnitCellForm(basis); // Niggli Second
  if (LDEBUG) {
    cerr << "WARNING: remeber to FixLattices for f2c,c2f and reciprocal and bring atoms in cell " << endl;
  }
  return basis;
}

// ******************************************************************************
// Fold atoms into cell
// ******************************************************************************
// Folds atoms into the cell

// xstructure::foldAtomsInCell()
// modify xstructure in-place
void xstructure::foldAtomsInCell(const xmatrix<double>& lattice_new, bool skew, double tol, bool check_min_dists) {
  // DX20210104
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " xstr(ORIG)=" << endl << (*this) << endl;
  }

  deque<_atom> atoms_new = ::foldAtomsInCell((*this), lattice_new, skew, tol, check_min_dists);
  // fold atoms in //DX20210118 - specify global namespace

  // update xstructure info
  (*this).lattice = lattice_new;
  // sort and update atom counts/order/types/basis/etc.
  std::stable_sort(atoms_new.begin(), atoms_new.end(), sortAtomsNames); // DX20210129
  (*this).ReplaceAtoms(atoms_new);
  (*this).BringInCell();
  (*this).FixLattices();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " xstr(NEW)=" << endl << (*this) << endl;
  }
}

deque<_atom> foldAtomsInCell(const xstructure& a, const xmatrix<double>& lattice_new, bool skew, double tol, bool check_min_dists) {
  // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 = added check_min_dists bool
  const bool LDEBUG = (false || XHOST.DEBUG);

  const double volume_original = std::abs(aurostd::det(a.lattice));
  const double volume_new = std::abs(aurostd::det(lattice_new));
  const bool fold_in_only = ((volume_new < volume_original) || (aurostd::isequal(volume_original, volume_new)));

  deque<_atom> atoms_orig = a.atoms; // need to make a copy for the pointer
  deque<_atom>* ptr_atoms = &atoms_orig;
  deque<_atom> atoms = *ptr_atoms;
  // DX20210129 - this need to be done before if-statement since we put atomic_grid inside if-statement
  if (!fold_in_only) {
    xstructure atomic_grid; // stays empty if not needed //DX+ME20210111 - added inside if-statement
    const double radius = RadiusSphereLattice(lattice_new);
    const xvector<int> dims = LatticeDimensionSphere(a.lattice, radius); // int dim=max(dims)+1; //dim=3;  //CO20190520
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " a.lattice=" << endl;
      cerr << a.lattice << endl;
      cerr << __AFLOW_FUNC__ << " lattice_new=" << endl;
      cerr << lattice_new << endl;
      cerr << __AFLOW_FUNC__ << " vol(a.lattice)=" << std::abs(aurostd::det(a.lattice)) << endl;
      cerr << __AFLOW_FUNC__ << " vol(lattice_new)=" << std::abs(aurostd::det(lattice_new)) << endl;
      cerr << __AFLOW_FUNC__ << " radius(a.lattice)=" << RadiusSphereLattice(a.lattice) << endl;
      cerr << __AFLOW_FUNC__ << " radius(lattice_new)=" << radius << endl;
      cerr << __AFLOW_FUNC__ << " dims=" << dims << endl;
    }
    //[CO20190520 - excessive, too large of an exploration radius]xmatrix<double> supercell; supercell(1,1)=dims(1); supercell(2,2)=dims(2); supercell(3,3)=dims(3);  //NO NEED, function ensures radius is encompassed //be safe and go +1 out
    // xmatrix<double> supercell; supercell(1,1)=dim; supercell(2,2)=dim; supercell(3,3)=dim;  //NO NEED, function ensures radius is encompassed //be safe and go +1 out
    // vector<int> sc2pcMap, pc2scMap; //dummy
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " building atomic grid with dims=[" << dims << "]" << endl;
    }
    atomic_grid = a;
    atomic_grid.clean(); // DX20191220 - uppercase to lowercase clean
    atomic_grid.GenerateGridAtoms(dims[1], dims[2], dims[3]); // much faster than supercell
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " atomic grid built" << endl;
    }
    ptr_atoms = &atomic_grid.grid_atoms; // CO20190808 - GenerateGridAtoms() populates grid_atoms, not atoms
    atoms = *ptr_atoms; // DX20210129 - set inside if-statement, otherwise, grid atoms goes out of scope
  }

  return foldAtomsInCell(atoms, a.lattice, lattice_new, skew, tol, check_min_dists);
  // DX20190619 = added check_min_dists bool
}

deque<_atom> foldAtomsInCell(const deque<_atom>& atoms,
                             const xmatrix<double>& lattice_orig,
                             const xmatrix<double>& lattice_new,
                             bool skew,
                             double tol,
                             bool check_min_dists) { // DX20190619 - added check_min_dists bool
  const bool LDEBUG = (false || XHOST.DEBUG);

  deque<_atom> atoms_in_cell;

  const xmatrix<double> f2c_new = trasp(lattice_new);
  const xmatrix<double> c2f_new = inverse(f2c_new);

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " f2c_new=" << endl;
    cerr << f2c_new << endl;
    cerr << __AFLOW_FUNC__ << " c2f_new=" << endl;
    cerr << c2f_new << endl;
  }

  _atom tmp;
  for (size_t j = 0; j < atoms.size(); j++) {
    //[CO20190520 - this case is not needed]if(atoms_in_cell.size() == 0) {
    //[CO20190520 - this case is not needed]  atoms_in_cell.push_back(atoms[j]);
    //[CO20190520 - this case is not needed]  atoms_in_cell.back().fpos = BringInCell(c2f_new * atoms[j].cpos);
    //[CO20190520 - this case is not needed]  atoms_in_cell.back().cpos = f2c_new * atoms_in_cell.back().fpos;
    //[CO20190520 - this case is not needed]  atoms_in_cell.back().ijk(1)=0; atoms_in_cell.back().ijk(2)=0; atoms_in_cell.back().ijk(3)=0;
    //[CO20190520 - this case is not needed]} else {  //[CO20200106 - close bracket for indenting]}
    // bool duplicate_atom = false;
    tmp.fpos = BringInCell(c2f_new * atoms[j].cpos);
    tmp.cpos = f2c_new * tmp.fpos;
    if (!SYM::MapAtom(atoms_in_cell, tmp, false, lattice_new, f2c_new, skew, tol)) {
      // DX20190619 - lattice_new and f2c_new as input
      atoms_in_cell.push_back(atoms[j]);
      atoms_in_cell.back().fpos = tmp.fpos;
      atoms_in_cell.back().cpos = tmp.cpos;
      atoms_in_cell.back().ijk(1) = 0;
      atoms_in_cell.back().ijk(2) = 0;
      atoms_in_cell.back().ijk(3) = 0;
    }
    //[CO20190520 - this case is not needed]}
  }

  if (check_min_dists) { // DX20190613
    const double min_dist_orig = SYM::minimumDistance(atoms);
    // lattice_orig //this does NOT work if we use GenerateGridAtoms (no longer periodic with lattice), so simply compare distances between atoms. NOTE: this is no longer the true minimumDistance(), which requires knowledge of the lattice vectors
    const double min_dist_new = SYM::minimumDistance(atoms_in_cell); // lattice_new
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lattice_orig=" << endl;
      cerr << lattice_orig << endl;
      cerr << __AFLOW_FUNC__ << " lattice_new=" << endl;
      cerr << lattice_new << endl;
      cerr << __AFLOW_FUNC__ << " atoms_orig=" << endl;
      for (size_t i = 0; i < atoms.size(); i++) {
        cerr << atoms[i] << endl;
      }
      cerr << __AFLOW_FUNC__ << " min_dist_orig=" << endl;
      cerr << min_dist_orig << endl;
      cerr << __AFLOW_FUNC__ << " min_dist_new=" << endl;
      cerr << min_dist_new << endl;
    }
    if (!aurostd::isequal(min_dist_orig, min_dist_new, 0.1)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Minimum distance changed, check that atoms are not rotated", _RUNTIME_ERROR_);
    }
  }

  // sort atoms //DX+CO20210119
  std::stable_sort(atoms_in_cell.begin(), atoms_in_cell.end(), sortAtomsNames);

  return atoms_in_cell;
}

// ***************************************************************************
// Function GetPrimitiveVASP()
// ***************************************************************************
// get reduced lattice in reciprocal space, converts to real and folds atoms into this lattice
// this should be the fastest lattice for VASP
// CO20180409 - refer to standard primitive instead, this function is simply a test, probably not optimal to standard primitive
xstructure GetPrimitiveVASP(const xstructure& a) {
  double tol = a.sym_eps;
  if (tol == AUROSTD_NAN) {
    tol = SYM::defaultTolerance(a);
  }
  return GetPrimitiveVASP(a, tol);
}

xstructure GetPrimitiveVASP(const xstructure& a, double tol) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  xstructure b(a);
  b.ReScale(1.0);
  xmatrix<double> klattice = ReciprocalLattice(b.lattice, b.scale); // get reciprocal lattice
  klattice = MinkowskiBasisReduction(klattice); // get min of lattice
  const xmatrix<double> lattice_new = ReciprocalLattice(klattice, b.scale); // get direct lattice
  // get skew - START
  if (b.dist_nn_min == AUROSTD_NAN) {
    b.dist_nn_min = SYM::minimumDistance(b.atoms, b.lattice);
  }
  const bool skew = SYM::isLatticeSkewed(lattice_new, b.dist_nn_min, tol);
  // get skew - STOP
  b.lattice = lattice_new; // b.f2c=trasp(b.lattice); b.c2f=inverse(b.f2c);  //set new lattice, f2c, c2f
  // b.atoms=foldAtomsInCell(b.atoms,b.c2f,b.f2c,skew,tol);                //fold atoms into new lattice
  //[CO20190520 - ReplaceAtoms()]b.atoms=foldAtomsInCell(a,b.lattice,skew,tol,true);                   //fold atoms into new lattice, don't bother folding out (slow)
  const deque<_atom> atoms = foldAtomsInCell(a, b.lattice, skew, tol); // CO20190520 - ReplaceAtoms()
  b.ReplaceAtoms(atoms); // CO20190520 - ReplaceAtoms()
  if ((false || LDEBUG) && !isequal(a.lattice, b.lattice, 1e-3)) {
    cerr << "-----------------------------------------------------------------------" << endl;
    cerr << "ORIG STRUCTURE " << endl;
    cerr << a;
    b.ReScale(1.0);
    b.ShiftOriginToAtom(0);
    b.BringInCell(); // fast clean for comparison
    cerr << "NEW STRUCTURE " << endl;
    cerr << b;
    cerr << "STRUCTURES IDENTICAL = " << (compare::structuresMatch(a, b, true) ? "true" : "false") << endl;
    cerr << "-----------------------------------------------------------------------" << endl;
  }
  return b;
}

// **************************************************************************
// Function BringInCell
// **************************************************************************
// these routines take atoms or structures and brings them
// to the unit cell (incell). There is some overloading for
// structures. SC Aug2007
// EDITED BY CO+DX to include a tolerance that converts from Cartesian
// space to direct space via covariant and contravariant transformations.
#define _EPS_roundoff_ 0.001
#define _incellcutoff_ (1.0 - _EPS_roundoff_)
// DX+CO, IF EPSILON IS PROVIDED, THEN WE ASSUME YOU WANT US TO MOVE THE ATOMS, OTHERWISE IT'S A HARD CUT OFF

// BringInCell() - ROBUST (DX+ME+CO20190905)
// this function brings components/xvectors/_atoms into a unit cell
// the upper bound and lower bound of the cell can be adjusted
// (e.g., standard unit cell : 0.0 to 1.0 ; unit cell centered on origin: -0.5 to 0.5)
// the tolerance can be tuned and "shifts" the bounds, favoring the lower bound
// (e.g., for bounds 0.0 to 1.0, we favor the origin, bringing values
// inside the cell if they are between "lower_bound-tolerance" and "upper_bound-tolerance")
// the AFLOW developers suggest a hard cutoff, e.g., _ZERO_TOL_ (DX+ME+CO)

// **************************************************************************
// BringInCellInPlace() (change value/object in place)

// -------------------------------------------------------------------
// double (change in place)
void BringInCellInPlace(double& component, double tolerance, double upper_bound, double lower_bound) {
  if (component == INFINITY || component != component || component == -INFINITY) {
    stringstream message; // Moving the stringstream outside the if-statement would add a lot to the run time (~1 sec).
    message << "Value of component is invalid: (+-) INF or NAN value (component=" << component << ").";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_); // DX20190905 - replaced cerr with throw
  }
  if (std::signbit(tolerance)) { // DX20191115
    stringstream message; // Moving the stringstream outside the if-statement would add a lot to the run time (~1 sec).
    message << "Sign of tolerance is negative (tolerance=" << tolerance << ").";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }
  while (component - upper_bound >= -tolerance) {
    component -= 1.0;
  }
  // note: non-symmetric, favors values closer to lower bound
  while (component - lower_bound < -tolerance) {
    component += 1.0;
  }
}

// -------------------------------------------------------------------
// xvector (change in place)
void BringInCellInPlace(xvector<double>& fpos, double tolerance, double upper_bound, double lower_bound) {
  for (int i = fpos.lrows; i <= fpos.urows; i++) {
    BringInCellInPlace(fpos[i], tolerance, upper_bound, lower_bound);
  }
}

// -------------------------------------------------------------------
// _atom (change in place, updates both fpos and pos)
void BringInCellInPlace(_atom& atom, const xmatrix<double>& lattice, double tolerance, double upper_bound,
                        double lower_bound) { // DX20190904
  BringInCellInPlaceFPOS(atom, tolerance, upper_bound, lower_bound); // update fpos first

  // update cpos
  atom.cpos = F2C(lattice, atom.fpos); // update cpos next
  atom.isincell = true;
}

// -------------------------------------------------------------------
// xstructure (change in place)
void BringInCellInPlace(xstructure& xstr, double tolerance, double upper_bound, double lower_bound) { // DX20190904
  for (size_t i = 0; i < xstr.atoms.size(); i++) {
    BringInCellInPlace(xstr.atoms[i], xstr.lattice, tolerance, upper_bound, lower_bound);
  }
}

// **************************************************************************
// BringInCell() (return new value/object)

// -------------------------------------------------------------------
// double (return new double)
double BringInCell(double component_in, double tolerance, double upper_bound, double lower_bound) {
  double component_out = component_in;
  if (component_out == INFINITY || component_out != component_out || component_out == -INFINITY) {
    stringstream message; // Moving the stringstream outside the if-statement would add a lot to the run time (~1 sec).
    message << "Value of component is invalid: (+-) INF or NAN value (component=" << component_out << ").";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_); // DX20190905 - replaced cerr with throw
  }
  if (std::signbit(tolerance)) { // DX20191115
    stringstream message; // Moving the stringstream outside the if-statement would add a lot to the run time (~1 sec).
    message << "Sign of tolerance is negative (tolerance=" << tolerance << ").";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }
  while (component_out - upper_bound >= -tolerance) {
    component_out -= 1.0;
  }
  // note: non-symmetric, favors values closer to lower bound
  while (component_out - lower_bound < -tolerance) {
    component_out += 1.0;
  }
  return component_out;
}

// -------------------------------------------------------------------
// xvector (return new xvector)
xvector<double> BringInCell(const xvector<double>& fpos_in, double tolerance, double upper_bound, double lower_bound) {
  // DX20190904
  const xvector<double> fpos_out = fpos_in;
  for (int i = fpos_out.lrows; i <= fpos_out.urows; i++) {
    BringInCellInPlace(fpos_out[i], tolerance, upper_bound, lower_bound);
  }
  return fpos_out;
}

// -------------------------------------------------------------------
// _atom (return new _atom, update fpos only)
_atom BringInCellFPOS(const _atom& atom_in, double tolerance, double upper_bound, double lower_bound) { // DX20190904
  _atom atom_out = atom_in;
  BringInCellInPlace(atom_out.fpos, tolerance, upper_bound, lower_bound);

  // update ijk
  for (int i = atom_out.fpos.lrows; i <= atom_out.fpos.urows; i++) {
    atom_out.ijk(i) = (int) (atom_out.fpos[i] - atom_in.fpos[i]);
  }
  return atom_out;
}

// -------------------------------------------------------------------
// xstructure (return xstructure)
xstructure BringInCell(const xstructure& xstr_in, double tolerance, double upper_bound, double lower_bound) {
  // DX20190904
  xstructure xstr_out = xstr_in;
  for (size_t i = 0; i < xstr_out.atoms.size(); i++) {
    BringInCellInPlace(xstr_out.atoms[i], xstr_out.lattice, tolerance, upper_bound, lower_bound);
  }
  return xstr_out;
}

// **************************************************************************
// BringInCellFPOS() (updates atom.fpos only)

// -------------------------------------------------------------------
// _atom (change in place, update fpos only)
void BringInCellInPlaceFPOS(_atom& atom, double tolerance, double upper_bound, double lower_bound) { // DX20190904
  const xvector<double> orig_fpos = atom.fpos; // DX - needed for ijk later
  BringInCellInPlace(atom.fpos, tolerance, upper_bound, lower_bound);

  // update ijk
  for (int i = atom.fpos.lrows; i <= atom.fpos.urows; i++) {
    atom.ijk(i) = (int) (atom.fpos[i] - orig_fpos[i]);
  }
}

// -------------------------------------------------------------------
// _atom (return new _atom, update fpos and cpos)
_atom BringInCell(const _atom& atom_in, const xmatrix<double>& lattice, double tolerance, double upper_bound,
                  double lower_bound) { // DX20190904
  _atom atom_out = BringInCellFPOS(atom_in, tolerance, upper_bound, lower_bound);

  // update cpos
  atom_out.cpos = F2C(lattice, atom_out.fpos);
  atom_out.isincell = true;

  return atom_out;
}

// **************************************************************************
// BringInCell() (method for xstructure)

// -------------------------------------------------------------------
// xstructure (xstructure method)
void xstructure::BringInCell(double tolerance, double upper_bound, double lower_bound) { // DX20190904
  for (size_t i = 0; i < atoms.size(); i++) {
    BringInCellInPlace(atoms[i], lattice, tolerance, upper_bound, lower_bound);
    // DX "::" to access outside of xstructure class
  }
}

// ***************************************************************************
// Function BringInCompact
// ***************************************************************************
// Make a structure where all the atoms are all the
// atoms are mapped through the unit and neighbors cells
// to minimize the shortest possible bond with an adjacent atom
// This option is very useful if you run big and complicate
// molecules where atoms exit of the unit cell and you have
// problems understanding where they are because visualization
// packages do not show bonds anymore ...
// Anyway, it is easier to test than to describe. (SC 6 Aug 04)
xstructure BringInCompact(const xstructure& a) {
  xstructure sstr = a;
  sstr.BringInCompact();
  return sstr;
}

void xstructure::BringInCompact() {
  double min_bond;
  double bond;
  // For direct coordinates
  BringInCell();
  xvector<double> adref1pos(3);
  const xvector<double> adtstpos(3);
  xvector<double> adtrgpos(3);
  const xvector<double> acref1pos(3);
  const xvector<double> actstpos(3);
  xvector<double> actrgpos(3);
  for (size_t i = 1; i < atoms.size(); i++) { // scan all the atoms to move except the first
    adtrgpos = atoms[i].fpos;
    actrgpos = atoms[i].cpos;
    min_bond = 1.0e6;
    for (uint ii = 0; ii < i; ii++) { // scan over all the reference atoms
      adref1pos = atoms[ii].fpos;
      for (int ic = 1; ic <= 3; ic++) { // scan over all the reference atoms
        acref1pos(ic) = 0.0; // scan over all the reference atoms
        for (int jc = 1; jc <= 3; jc++) {
          acref1pos(ic) = acref1pos(ic) + adref1pos(jc) * lattice(jc, ic);
        }
        // scan over all the reference atoms
      } // scan over all the reference atoms

      for (int i1 = -1; i1 <= 1; i1++) { // roll over first neighbor cells
        for (int j1 = -1; j1 <= 1; j1++) { // roll over first neighbor cells
          for (int k1 = -1; k1 <= 1; k1++) { // roll over first neighbor cells
            adtstpos(1) = atoms[i].fpos(1) + i1;
            adtstpos(2) = atoms[i].fpos(2) + j1;
            adtstpos(3) = atoms[i].fpos(3) + k1; // test the atom
            for (int ic = 1; ic <= 3; ic++) { // test the atom
              actstpos(ic) = 0.0; // test the atom
              for (int jc = 1; jc <= 3; jc++) { // test the atom
                actstpos(ic) = actstpos(ic) + adtstpos(jc) * lattice(jc, ic); // test the atom
              }
            } // test the atom
            bond = aurostd::modulus(actstpos - acref1pos); // test the bond distance
            // if((bond<min_bond)                  // if it is OK then DO IT !
            if ((bond < 1.03 * min_bond && std::abs((int) i - (int) ii) < 10) || (bond < 0.98 * min_bond)) {
              // CO20200106 - patching for auto-indenting
              //  if it is OK then DO IT !
              min_bond = bond; // update
              adtrgpos = adtstpos;
              actrgpos = actstpos;
            }
          }
        }
      }
    }
    // now the TRG atom has all the information about the bst position of the atom to be moved.
    // than EAT IT !
    for (int j = 1; j <= 3; j++) {
      atoms[i].fpos(j) = adtrgpos(j);
      atoms[i].cpos(j) = actrgpos(j);
    }
  }
  // For Cartesian coordinates (get from direct coords);
  for (size_t iat = 0; iat < atoms.size(); iat++) {
    atoms[iat].cpos = F2C(lattice, atoms[iat].fpos);
  }
}

// ***************************************************************************
// Function BringInWignerSeitz
// ***************************************************************************
// Make a structure where all the atoms are
// mapped to their images in the Wigner-Seitz cell.(SC 10Jan04)
void xstructure::BringInWignerSeitz() {
  xvector<double> rat(3);
  xvector<double> a1(3);
  xvector<double> a2(3);
  xvector<double> a3(3);
  xvector<double> a12(3);
  xvector<double> a31(3);
  xvector<double> a23(3);
  xvector<double> a123(3);
  double na1;
  double na2;
  double na3;
  double na12;
  double na31;
  double na23;
  double na123;
  a1 = lattice(1);
  na1 = aurostd::modulus(a1);
  a2 = lattice(2);
  na2 = aurostd::modulus(a2);
  a3 = lattice(3);
  na3 = aurostd::modulus(a3);
  a12 = lattice(1) + lattice(2);
  na12 = aurostd::modulus(a12);
  a31 = lattice(1) + lattice(3);
  na31 = aurostd::modulus(a31);
  a23 = lattice(2) + lattice(3);
  na23 = aurostd::modulus(a23);
  a123 = lattice(1) + lattice(2) + lattice(3);
  na123 = aurostd::modulus(a123);
  double proj_a1;
  double proj_a2;
  double proj_a3;
  double proj_a12;
  double proj_a31;
  double proj_a23;
  double proj_a123;

  BringInCell();

  for (size_t iat = 0; iat < atoms.size(); iat++) {
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          rat = atoms[iat].cpos + ((double) i) * a1 + ((double) j) * a2 + ((double) k) * a3;
          proj_a1 = scalar_product(rat, a1) / na1 / na1;
          proj_a2 = scalar_product(rat, a2) / na2 / na2;
          proj_a3 = scalar_product(rat, a3) / na3 / na3;
          proj_a12 = scalar_product(rat, a12) / na12 / na12;
          proj_a31 = scalar_product(rat, a31) / na31 / na31;
          proj_a23 = scalar_product(rat, a23) / na23 / na23;
          proj_a123 = scalar_product(rat, a123) / na123 / na123;
          if ((proj_a1 > -0.5 && proj_a1 <= 0.5) && (proj_a2 > -0.5 && proj_a2 <= 0.5) && (proj_a3 > -0.5 && proj_a3 <= 0.5) && (proj_a12 > -0.5 && proj_a12 <= 0.5) && (proj_a31 > -0.5 && proj_a31 <= 0.5) &&
              (proj_a23 > -0.5 && proj_a23 <= 0.5) && (proj_a123 > -0.5 && proj_a123 <= 0.5)) {
            atoms[iat].cpos(1) = rat(1);
            atoms[iat].cpos(2) = rat(2);
            atoms[iat].cpos(3) = rat(3);
            i = 10;
            j = 10;
            k = 10;
          }
        }
      }
    }
    atoms[iat].fpos = C2F(lattice, atoms[iat].cpos);
  }
}

xstructure BringInWignerSeitz(const xstructure& a) {
  xstructure str = a;
  str.BringInWignerSeitz();
  return str;
}

// **************************************************************************
// Function GetPrimitive
// Function IsTranslationFVector& IsTranslationCVector
// ***************************************************************************
// This function returns 1 if the given vector is a translation
// vector of the given structure and 0 otherwise. Input translation
// vector is expectd to be in fractional coordinates.
// written by Stefano Curtarolo (superfast edition)
// #define IsTranslationFVector IsTranslationFVectorFAST
// #define IsTranslationFVector IsTranslationFVectorORIGINAL
// #define IsTranslationFVector IsTranslationFVectorFAST_2011
#define IsTranslationFVector IsTranslationFVectorORIGINAL_2011

// isTranslationVector() //DX20210316 - moved from XtalFinder
// faster than subsequent variants below and more robust
bool isTranslationVector(const xstructure& xstr, const xvector<double>& vec, double tolerance, bool is_frac) {
  // Check if input vector is a translation vector (i.e. preserves periodicity)
  // tolerance: default is half an Angstrom
  // (tolerance example: need at least tol=0.1 for As1_ICSD_158474 == As1_ICSD_162840 via XtalFinder)

  if (tolerance < _ZERO_TOL_) {
    stringstream message;
    message << "Zero tolerance: " << tolerance;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }

  xvector<double> cvec;
  xvector<double> fvec;
  if (is_frac) {
    fvec = vec;
    cvec = xstr.f2c * vec;
  } else {
    fvec = xstr.c2f * vec;
    cvec = vec;
  }

  const size_t natoms = xstr.atoms.size();
  const bool skew = false;
  uint count = 0;

  // ---------------------------------------------------------------------------
  // check if applying the translation maps to another atom
  // use MapAtom to match type/name/spin/occupation/etc.
  for (uint d = 0; d < natoms; d++) {
    _atom tmp_atom = xstr.atoms[d];
    tmp_atom.cpos += cvec;
    tmp_atom.fpos += fvec;
    if (SYM::MapAtom(xstr.atoms, tmp_atom, true, xstr.lattice, xstr.f2c, skew, tolerance)) {
      count++;
    }
    // match not found, violates periodicity, return immediately
    else {
      return false;
    }
  }
  return (count == natoms);
}

bool IsTranslationFVectorFAST(const xstructure& a, const xvector<double>& ftvec) {
  stringstream message;
  if (a.equiv_fpos_epsilon < 1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
  const double tolerance = a.equiv_fpos_epsilon;
  if (aurostd::modulus(ftvec) <= tolerance) {
    return true;
  }

  uint CntGoodTrans = 0;
  const xvector<double> ftpos(3);
  xvector<double> diff(3);
  bool found_atom = false;
  bool found_tvec = true;
  const xvector<int> types(0, a.atoms.size()); // mirror to enhance speed !
  for (size_t iat = 0; iat < a.atoms.size(); iat++) { // mirror to enhance speed !
    types(iat) = a.atoms[iat].type; // mirror to enhance speed !
  }
  for (size_t iat = 0; iat < a.atoms.size() && found_tvec; iat++) {
    // Get translated position and shift back into unit cell.
    //  ftpos=a.atoms.at(iat).fpos+ftvec;
    // ftpos=BringInCell(ftpos,tolerance);  // no because it will be needed later...
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    found_atom = false;
    for (size_t jat = 0; jat < a.atoms.size() && !found_atom && found_tvec; jat++) {
      if (types[iat] == types[jat]) {
        //        diff=ftpos-a.atoms[jat].fpos;
        diff = BringInCell(a.atoms[iat].fpos - a.atoms[jat].fpos + ftvec, tolerance);
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if (aurostd::modulus(diff) <= tolerance) {
          CntGoodTrans++;
          found_atom = true;
        }
      }
    } // For jat
    if (found_atom == false) {
      found_tvec = false; // at least one atom does not have a copy...
    }
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  // cout << "GOODs=" << CntGoodTrans << endl;
  if (CntGoodTrans == a.atoms.size()) {
    return true;
  }
  return false;
}

bool IsTranslationFVectorORIGINAL(const xstructure& a, const xvector<double>& ftvec) {
  stringstream message;
  if (a.equiv_fpos_epsilon < 1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
  const double tolerance = a.equiv_fpos_epsilon;
  if (aurostd::modulus(ftvec) <= tolerance) {
    return true;
  }
  uint CntGoodTrans = 0;
  const xvector<double> ftpos(3);
  xvector<double> diff(3);
  const xvector<int> types(0, a.atoms.size()); // mirror to enhance speed !
  for (size_t iat = 0; iat < a.atoms.size(); iat++) { // mirror to enhance speed !
    types(iat) = a.atoms[iat].type; // mirror to enhance speed !
  }
  for (size_t iat = 0; iat < a.atoms.size(); iat++) {
    // Get translated position and shift back into unit cell.
    // ftpos=a.atoms[iat].fpos+ftvec;
    // ftpos=BringInCell(ftpos,tolerance); // no because it will be needed later...
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    for (size_t jat = 0; jat < a.atoms.size(); jat++) {
      if (types[iat] == types[jat]) {
        //        diff=ftpos-a.atoms[jat].fpos;
        diff = BringInCell(a.atoms[iat].fpos - a.atoms[jat].fpos + ftvec, tolerance);
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if (aurostd::modulus(diff) <= tolerance) {
          CntGoodTrans++;
        }
      }
    } // For jat
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  //  cout << "GOODs=" << CntGoodTrans << endl;
  if (CntGoodTrans == a.atoms.size()) {
    return true;
  }
  return false;
}

bool IsTranslationFVectorFAST_2011(const xstructure& a, const xvector<double>& ftvec) {
  stringstream message;
  if (a.equiv_fpos_epsilon < 1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
  const double tolerance = a.equiv_fpos_epsilon;
  if (aurostd::modulus(ftvec) <= tolerance) {
    return true;
  }
  uint CntGoodTrans = 0;
  xvector<double> ftpos(3);
  xvector<double> diff(3);
  bool found_atom = false;
  bool found_tvec = true;
  const xvector<int> types(0, a.atoms.size()); // mirror to enhance speed !
  for (size_t iat = 0; iat < a.atoms.size(); iat++) { // mirror to enhance speed !
    types(iat) = a.atoms[iat].type; // mirror to enhance speed !
  }
  for (size_t iat = 0; iat < a.atoms.size() && found_tvec; iat++) {
    // Get translated position and shift back into unit cell.
    ftpos = a.atoms[iat].fpos + ftvec;
    ftpos = BringInCell(ftpos);
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    found_atom = false;
    for (size_t jat = 0; jat < a.atoms.size() && !found_atom && found_tvec; jat++) {
      if (types[iat] == types[jat]) {
        diff = ftpos - a.atoms[jat].fpos;
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if (aurostd::modulus(diff) < tolerance) {
          CntGoodTrans++;
          found_atom = true;
        }
      }
    } // For jat
    if (found_atom == false) {
      found_tvec = false; // at least one atom does not have a copy...
    }
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  if (CntGoodTrans == a.atoms.size()) {
    return true;
  }
  return false;
}

bool IsTranslationFVectorORIGINAL_2011(const xstructure& a, const xvector<double>& ftvec) {
  stringstream message;
  if (a.equiv_fpos_epsilon < 1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
  const double tolerance = a.equiv_fpos_epsilon;
  if (aurostd::modulus(ftvec) <= tolerance) {
    return true;
  }

  uint CntGoodTrans = 0;
  xvector<double> ftpos(3);
  xvector<double> diff(3);
  const xvector<int> types(0, a.atoms.size()); // mirror to enhance speed !
  for (size_t iat = 0; iat < a.atoms.size(); iat++) { // mirror to enhance speed !
    types(iat) = a.atoms[iat].type; // mirror to enhance speed !
  }
  for (size_t iat = 0; iat < a.atoms.size(); iat++) {
    // Get translated position and shift back into unit cell.
    ftpos = a.atoms[iat].fpos + ftvec;
    ftpos = BringInCell(ftpos);
    // Find closest atom in unit cell to translated position (should only be one atom less than tolerance).
    for (size_t jat = 0; jat < a.atoms.size(); jat++) {
      if (types[iat] == types[jat]) {
        diff = ftpos - a.atoms[jat].fpos;
        diff = SYM::minimizeDistanceFractionalMethod(diff); // DX20190613
        // If the translated atom maps onto another and its type is the same as the atom mapped onto then increment.
        if (aurostd::modulus(diff) < tolerance) {
          CntGoodTrans++;
        }
      }
    } // For jat
  } // For iat
  // If we counted a good translation for every atoms then we have a translation vector.
  if (CntGoodTrans == a.atoms.size()) {
    return true;
  }
  return false;
}

bool IsTranslationCVector(const xstructure& a, const xvector<double>& ctvec) {
  stringstream message;
  if (a.equiv_fpos_epsilon < 1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
  /* Input translation vector is expectd to be in cartesian coordinates. */
  return IsTranslationFVector(a, C2F(a.lattice, ctvec));
}

// **************************************************************************
// Function GetPrimitive
// **************************************************************************
// This funtion returns a structure object where everything
// has been transformed into a primitive lattice.
// Algorithm
// Construct a list of all candidate primitive lattice vectors.
// This list consists of the present lattice vectors and all basis
// vectors where translation by that basis vector leaves all the
// basis atom types unchanged. Loop over all possible triads in this list.
// For each triad find the volume. Store a list of all triads with the
// minimum volume. Then select the triad for which the first triad vector
// has the maximal projection onto the first lattice vector.  If there
// are more than one of these do the same for the second lattice vector,
// and then the third if needed. Then just pick one if there are still
// multiple candidates.

xstructure GetPrimitiveSINGLE(const xstructure& a, double tolerance);
xstructure GetPrimitiveMULTITHREAD(const xstructure& a, double tolerance);
xstructure GetPrimitive(const xstructure& a);

xstructure GetPrimitive(const xstructure& _a, double tolerance) {
  // return GetPrimitiveMULTITHREAD(a);

  xstructure b;
  b = GetPrimitive_20210322(_a, tolerance); // DX20210406 - new/fast variant

  //[DX TESTING] xstructure c;
  //[DX TESTING] c=GetPrimitiveMULTITHREAD(a,tolerance);
  //[DX TESTING] cerr << "STRUCTURES MATCH (ORIG vs OLD): " << compare::structuresMatch(a,b,false) << " (ORIG vs NEW): " << compare::structuresMatch(a,c,false) << endl;
  //[DX TESTING] cerr << "VOLUMES: " << b.atoms.size()/c.atoms.size() << " (mod=" << b.atoms.size()%c.atoms.size() << ")" << endl;

  // OLD VERSION - the version above is fast and does not need to be multithreaded

  return b;
}

xstructure GetPrimitive(const xstructure& a) {
  double tolerance = AUROSTD_MAX_DOUBLE;
  if (a.sym_eps != AUROSTD_NAN) {
    tolerance = a.sym_eps;
  } else {
    tolerance = SYM::defaultTolerance(a);
  }
  return GetPrimitive(a, tolerance);
}

// **************************************************************************
// GetPrimitive() //DX20210406
// **************************************************************************
// This version is faster than the previous AFLOW variants
// Speed ups:
//   - returns immediately if the number of atoms/types indicate the cell
//     cannot be reduced (i.e., only one atom or one atom of a given type)
//   - optimized lattice vector search (search over only the least-frequent
//     atom type)
//   - optimized lattice search (volume checks, mulitplicity checks, etc.)
//   - uses transfomation method to convert between original and primitive
//     cell (as opposed to removing atoms)
// It also uses sym_eps to reduce the cell (previous methods did not)
xstructure GetPrimitive_20210322(const xstructure& a, double eps) { // DX20210406
  xstructure xstr = a;
  xstr.GetPrimitive_20210322(eps);
  return xstr;
}

void xstructure::GetPrimitive_20210322(double eps) { // DX20210406

  const bool LDEBUG = (false || XHOST.DEBUG);

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN " << endl;
  }

  const size_t natoms_orig = (*this).atoms.size();
  // ---------------------------------------------------------------------------
  // if only one atom or atom type in the unit cell, the structure is already
  // primitivized
  if (natoms_orig == 1 || aurostd::min((*this).num_each_type) == 1) {
    return;
  }

  (*this).ReScale(1.0);
  (*this).FixLattices(); // DX20210407 - since we use c2f/f2c, update for safety

  // ---------------------------------------------------------------------------
  // set tolerance
  double tolerance = eps;
  if ((*this).dist_nn_min == AUROSTD_NAN) {
    (*this).MinDist();
  }

  if (tolerance == AUROSTD_MAX_DOUBLE) {
    if ((*this).sym_eps != AUROSTD_NAN) {
      tolerance = (*this).sym_eps;
    } else {
      tolerance = SYM::defaultTolerance((*this));
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [1] Tolerance = " << tolerance << endl;
  }

  const double volume_orig = (*this).Volume();

  // ---------------------------------------------------------------------------
  // get least frequent atom type and the corresponding set of atoms
  // to search for possible lattice vectors (minimal set of atoms perserving
  // periodicity)
  const uint atom_type_min = getLeastFrequentAtomTypes((*this))[0];
  // normally a vector, grabbing only first one (there will always be one)
  vector<uint> vindices_atoms_min = getAtomIndicesByType((*this), atom_type_min);
  const size_t natoms_min = vindices_atoms_min.size();

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [2] Subset of atoms to find lattice vectors: " << natoms_min << endl;
  }

  // generate list of vectors
  vector<xvector<double>> candidate_lattice_vector;
  candidate_lattice_vector.push_back((*this).lattice(1)); // lattice is made of good vectors
  candidate_lattice_vector.push_back((*this).lattice(2)); // lattice is made of good vectors
  candidate_lattice_vector.push_back((*this).lattice(3)); // lattice is made of good vectors

  // ---------------------------------------------------------------------------
  // get all lattice vectors
  // only need to check difference between 0th and ith atom coordinates for the
  // subset of atom indices
  vector<xvector<double>> diff_vectors;
  for (uint i = 1; i < natoms_min; i++) {
    diff_vectors.push_back(::BringInCell((*this).atoms[vindices_atoms_min[i]].fpos - (*this).atoms[vindices_atoms_min[0]].fpos));
    // need to get BringInCell from outside xstructure scope
  }
  // Translate by difference vectors and check if equivalent
  const bool is_frac = true;
  for (size_t d = 0; d < diff_vectors.size(); d++) {
    if (isTranslationVector((*this), diff_vectors[d], tolerance, is_frac)) {
      candidate_lattice_vector.push_back((*this).f2c * diff_vectors[d]);
    }
  }

  const double nlattice_vectors = candidate_lattice_vector.size();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " [3] number of lattice vectors=" << nlattice_vectors << endl;
  }

  if (LDEBUG) {
    for (uint i = 0; i < nlattice_vectors; i++) {
      cerr << __AFLOW_FUNC__ << i << ": mod=" << aurostd::modulus(candidate_lattice_vector[i]) << " vec=" << candidate_lattice_vector[i] << endl;
    }
  }

  // ---------------------------------------------------------------------------
  // if only the three original lattice vectors were found, then the cell
  // is already primitivized
  if (nlattice_vectors == 3) {
    return;
  }

  const xmatrix<double> lattice_tmp(3, 3);
  xmatrix<double> plattice(3, 3);
  xmatrix<double> olattice(3, 3);
  olattice = (*this).lattice; // the lattice is always a good lattice
  double volume_min = aurostd::det((*this).lattice);
  double volume_tmp = AUROSTD_MAX_DOUBLE;
  const double amin = (*this).dist_nn_min;
  const double tol_vol = tolerance * amin * amin * amin; // scale tolerance for volume
  double atom_number_ratio = 1; // track atom ratio during reduction
  const double integer_tol = 0.05; // use somewhat large tol to find lattices
  double moduli_sum = AUROSTD_MAX_DOUBLE;
  double moduli_sum_tmp = AUROSTD_MAX_DOUBLE;
  // ensures shortest sum of lattice vectors (chooses minimum Minkowski lattice, some noise can get introduced)

  // ---------------------------------------------------------------------------
  // check possible lattices, upper triangular for-loop only
  // criteria for volume of new lattice:
  //  1) must not be zero
  //  2) must be smaller than original lattice volume
  //  3) must be an integer multiple of original lattice
  //  4) must reduce atom count consistent with integer multiple
  //  5) must not remove entire atom types
  //  6) must form shortest vectors
  // Minkowski/Niggli will fix left-handed (negative determinants)
  for (uint iu = 0; iu < nlattice_vectors; iu++) {
    for (uint i = 1; i <= 3; i++) {
      lattice_tmp[1][i] = candidate_lattice_vector.at(iu)[i];
    }
    for (uint iv = iu + 1; iv < nlattice_vectors; iv++) {
      for (uint i = 1; i <= 3; i++) {
        lattice_tmp[2][i] = candidate_lattice_vector.at(iv)[i];
      }
      for (uint iw = iv + 1; iw < nlattice_vectors; iw++) {
        for (uint i = 1; i <= 3; i++) {
          lattice_tmp[3][i] = candidate_lattice_vector.at(iw)[i];
        }
        volume_tmp = std::abs(lattice_tmp[1][1] * lattice_tmp[2][2] * lattice_tmp[3][3] + lattice_tmp[1][2] * lattice_tmp[2][3] * lattice_tmp[3][1] + // FAST
                              lattice_tmp[1][3] * lattice_tmp[2][1] * lattice_tmp[3][2] - lattice_tmp[1][3] * lattice_tmp[2][2] * lattice_tmp[3][1] - // FAST
                              lattice_tmp[1][2] * lattice_tmp[2][1] * lattice_tmp[3][3] - lattice_tmp[1][1] * lattice_tmp[2][3] * lattice_tmp[3][2]); // FAST
        if (aurostd::abs(volume_tmp) > tol_vol) {
          plattice = ::MinkowskiBasisReduction(lattice_tmp);
          // Minkowski first, "::" is needed to access outside xstructure scope
          plattice = ::NiggliUnitCellForm(plattice); // Niggli Second, "::" is needed to access outside xstructure scope
          volume_tmp = (plattice[1][1] * plattice[2][2] * plattice[3][3] + plattice[1][2] * plattice[2][3] * plattice[3][1] + // FAST
                        plattice[1][3] * plattice[2][1] * plattice[3][2] - plattice[1][3] * plattice[2][2] * plattice[3][1] -
                        // FAST
                        plattice[1][2] * plattice[2][1] * plattice[3][3] - plattice[1][1] * plattice[2][3] * plattice[3][2]);
          // FAST
          atom_number_ratio = (double) natoms_min / aurostd::nint(volume_orig / volume_tmp);
          moduli_sum_tmp = aurostd::modulus(plattice(1)) + aurostd::modulus(plattice(2)) + aurostd::modulus(plattice(3));
          if (volume_tmp > tol_vol && // 1) tol > 0
              volume_tmp < volume_orig && // 2) new_tol < orig_tol
              aurostd::isinteger(volume_orig / volume_tmp, integer_tol) &&
              // 3) new and orig volumes related by integer multiple
              aurostd::isinteger(atom_number_ratio, 0.05) && // 4) number of atoms is consistent with integer multiple
              (atom_number_ratio - 1.0) > -_ZERO_TOL_ && // 5) didn't lose any atoms
              moduli_sum_tmp < moduli_sum) {
            // 6) ensure shortest sum of vectors (Minkowski should ensure this, but noise can be introduced)
            moduli_sum = moduli_sum_tmp;
            olattice = plattice;
            volume_min = (olattice[1][1] * olattice[2][2] * olattice[3][3] + olattice[1][2] * olattice[2][3] * olattice[3][1] + // FAST
                          olattice[1][3] * olattice[2][1] * olattice[3][2] - olattice[1][3] * olattice[2][2] * olattice[3][1] -
                          // FAST
                          olattice[1][2] * olattice[2][1] * olattice[3][3] - olattice[1][1] * olattice[2][3] * olattice[3][2]);
            // FAST
          }
        }
      }
    }
  }
  plattice = olattice;

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " reduced lattice=" << plattice << endl;
  }

  // ---------------------------------------------------------------------------
  // if the lattice remains the same, do not change or update the atoms
  if (aurostd::isequal(volume_min, volume_orig) || aurostd::isequal(plattice, (*this).lattice)) {
    return;
  }

  // ---------------------------------------------------------------------------
  // DX20210316 - used transformation method (more efficient)
  // ---------------------------------------------------------------------------
  // update xstructure
  xstructure prim = (*this);
  // prim.FixLattices();
  const xmatrix<double> transformation_matrix = GetBasisTransformation((*this).lattice, plattice);
  const xmatrix<double> rotation_matrix = aurostd::eye<double>(3, 3);
  // try to primitivize, it may fail, so return original structure
  try {
    // this checks volumes, number of atoms, etc. internally
    prim.TransformStructure(transformation_matrix, rotation_matrix);
  } catch (aurostd::xerror& err) { // CO20230220 - patching names vs. types bug
    if (aurostd::substring2bool(err.whereFileName(), "aflow_xatom.cpp") && aurostd::substring2bool(err.whereFunction(), "xstructure::ReplaceAtoms():") &&
        aurostd::substring2bool(err.what(), "atoms and species were incorrectly (re)sorted")) {
      throw err;
    }
    return;
  }

  // no messed up volume
  // this is checked in TransformStructure, but kept as a safety
  // DX20210623 - the check has been improved:
  // first, check if the reduction factor (inverse of fraction) is an integer
  // then, check if that factor is consistent with the number of atoms
  // this method is less sensitive to the tolerance threshold
  const double fraction = aurostd::det(prim.lattice) / aurostd::det((*this).lattice);
  const double reduction_factor = 1.0 / fraction;
  if (!aurostd::isinteger(reduction_factor, 0.1)) {
    stringstream message;
    message << "The original volume is not an integer multiple of the new volume: " << reduction_factor;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }
  const uint reduction_factor_integer = (uint) round(reduction_factor);
  if (std::abs(natoms_orig - (double) reduction_factor_integer * prim.atoms.size()) > 0.1) {
    stringstream message;
    message << "ERROR   " << __AFLOW_FUNC__ << endl;
    message << "        supercell has the wrong number of atoms" << endl;
    message << "        volume original    = " << (*this).Volume() << endl;
    message << "        volume prim        = " << prim.Volume() << endl;
    message << "        a.scale            = " << (*this).scale << endl;
    message << "        b.scale            = " << prim.scale << endl;
    message << "        a.atoms.size()     = " << (*this).atoms.size() << endl;
    message << "        b.atoms.size()     = " << prim.atoms.size() << endl;
    message << "        fraction           = " << fraction << endl;
    message << "        reduction_factor   = " << reduction_factor << endl; // DX20210623
    message << "        supercell atoms    = " << reduction_factor * prim.atoms.size() << endl; // DX20210623
    message << prim << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }
  // everything ok
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END [ok]=" << fraction << endl;
  }

  // set primitive representation
  (*this) = prim;
}

xstructure GetPrimitiveSINGLE(const xstructure& _a, double tolerance) { // APRIL 2009JUNE 2012 added tolerance
  const bool LDEBUG = (false || XHOST.DEBUG);
  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);
  xstructure a(_a);
  xstructure sstr = a;
  sstr.SetVolume(sstr.atoms.size());
  sstr = ReScale(sstr, 1.0);
  sstr = BringInCell(sstr);
  if (tolerance <= 0.0) {
    a.equiv_fpos_epsilon = _EQUIV_FPOS_EPS_;
  } else {
    a.equiv_fpos_epsilon = tolerance;
  }
  const double sstr_volume = sstr.Volume();

  const _aflags aflags;

  xmatrix<double> plattice(3, 3);
  xmatrix<double> olattice(3, 3);
  xvector<double> fdisp(3);
  xvector<double> cdisp(3);
  std::vector<xvector<double>> candidate_lattice_vector;

  int specie_min = 100000;
  int ispecie_min = 0;
  for (size_t ispecie = 0; ispecie < sstr.num_each_type.size(); ispecie++) {
    if (sstr.num_each_type[ispecie] < specie_min) {
      specie_min = sstr.num_each_type[ispecie];
      ispecie_min = ispecie;
    }
  }

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1)); // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2)); // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3)); // lattice is made of good vectors
  for (size_t iat1 = 0; iat1 < sstr.atoms.size(); iat1++) {
    if (sstr.atoms[iat1].type == ispecie_min) {
      for (size_t iat2 = 0; iat2 < sstr.atoms.size(); iat2++) {
        if (sstr.atoms[iat2].type == ispecie_min) {
          fdisp = sstr.atoms[iat2].fpos - sstr.atoms[iat1].fpos;
          cdisp = sstr.atoms[iat2].cpos - sstr.atoms[iat1].cpos;
          if (aurostd::modulus(fdisp) > 0.01 && aurostd::modulus(cdisp) > 0.01) {
            if (IsTranslationFVector(sstr, fdisp)) {
              candidate_lattice_vector.push_back(cdisp);
            }
          }
        }
      }
    }
  }
  if (LDEBUG) {
    cout << "DEBUG" << candidate_lattice_vector.size() << endl;
  }
  int cnt = 0;
  olattice = sstr.lattice; // the lattice is always a good lattice
  // now generate triplets
  for (size_t iu = 0; iu < candidate_lattice_vector.size(); iu++) {
    for (uint i = 1; i <= 3; i++) {
      plattice(1, i) = candidate_lattice_vector.at(iu)(i);
    }
    for (size_t iv = 0; iv < candidate_lattice_vector.size() && iv != iu; iv++) {
      for (uint i = 1; i <= 3; i++) {
        plattice(2, i) = candidate_lattice_vector.at(iv)(i);
      }
      for (size_t iw = 0; iw < candidate_lattice_vector.size() && iw != iv && iw != iu; iw++) {
        for (uint i = 1; i <= 3; i++) {
          plattice(3, i) = candidate_lattice_vector.at(iw)(i);
        }
        if (det(plattice) > 0.999 && det(plattice) < sstr_volume) {
          // no coplanar and contain at least 1 atom and smaller than the original cell
          if (aurostd::isinteger(sstr_volume / det(plattice))) { // integer ratio of volumes
            if (det(plattice) < det(olattice)) { // better than before
              if (LDEBUG) {
                cout << XPID << "DEBUG" << iu << "," << iv << "," << iw << " " << sstr_volume << " " << det(plattice) << " " << sstr_volume / det(plattice) << endl;
              }
              if (isdifferent(plattice, olattice, 0.0001)) {
                plattice = MinkowskiBasisReduction(plattice); // Minkowski first
                plattice = NiggliUnitCellForm(plattice); // Niggli Second
                if (isdifferent(plattice, olattice, 0.0001)) {
                  olattice = plattice;
                  cnt++;
                }
              }
            }
          }
        }
      }
    }
  }
  plattice = olattice;
  // done

  xstructure b = sstr;
  b.lattice = plattice; // b.lattice=roundoff(b.lattice,_EPS_FPOS_EQUAL_);
  b.FixLattices();
  b.write_lattice_flag = false;
  b.write_klattice_flag = false;
  b.write_DEBUG_flag = false;
  // plug them all
  for (size_t iat = 0; iat < sstr.atoms.size(); iat++) {
    b.atoms.at(iat).fpos = BringInCell(C2F(b.lattice, b.atoms.at(iat).cpos));
    b.atoms.at(iat).cpos = F2C(b.lattice, b.atoms.at(iat).fpos);
  }
  // now remove them
  b.RemoveFractionalCopies();
  // rescale back to original scale.
  b.SetVolume(Volume(a) * b.atoms.size() / a.atoms.size());
  b = ReScale(b, a.scale);
  //  // fix it up with the new Minkowsky and Niggli reductions // CANT DO AUTOMATICALLY
  //  b=LatticeReduction(b);
  // Put everything in new primitive cell.
  b = BringInCell(b);
  // check !
  const double fraction = Volume(a) / Volume(b);
  if (std::abs(b.atoms.size() * fraction - a.atoms.size()) > 0.1) {
    stringstream message;
    message << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    message << "        supercell has the wrong number of atoms" << endl;
    message << "        volume original    = " << Volume(a) << endl;
    message << "        volume prim        = " << Volume(b) << endl;
    message << "        a.scale            = " << a.scale << endl;
    message << "        b.scale            = " << b.scale << endl;
    message << "        a.atoms.size()     = " << a.atoms.size() << endl;
    message << "        b.atoms.size()     = " << b.atoms.size() << endl;
    message << "        fraction           = " << fraction << endl;
    message << "        supercell atoms    = " << fraction * b.atoms.size() << endl;
    message << b << endl;
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::GetPrimitiveSINGLE()", message, _RUNTIME_ERROR_);
  }
  // everything ok
  b.ClearSymmetry(); // CO20181226 - new structure, symmetry not calculated
  b.primitive_calculated = true; // DX20201007
  return b;
}

xstructure GetPrimitive1(const xstructure& a) { // MARCH 2009
  stringstream message;
  if (a.equiv_fpos_epsilon < 1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
  const double tolerance = a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);
  xstructure sstr = a;
  sstr.SetVolume(sstr.atoms.size());
  sstr = ReScale(sstr, 1.0);
  sstr = BringInCell(sstr);
  const double sstr_volume = sstr.Volume();
  const _aflags aflags;
  // identify the minimum set of atoms

  xmatrix<double> plattice(3, 3);
  xmatrix<double> olattice(3, 3);
  xvector<double> fdisp(3);
  xvector<double> cdisp(3);
  std::vector<xvector<double>> candidate_lattice_vector;

  int specie_min = 100000;
  int ispecie_min = 0;
  for (size_t ispecie = 0; ispecie < sstr.num_each_type.size(); ispecie++) {
    if (sstr.num_each_type[ispecie] < specie_min) {
      specie_min = sstr.num_each_type[ispecie];
      ispecie_min = ispecie;
    }
  }

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1)); // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2)); // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3)); // lattice is made of good vectors
  olattice = sstr.lattice; // the lattice is always a good lattice
  for (size_t iat1 = 0; iat1 < sstr.atoms.size(); iat1++) {
    for (size_t iat2 = 0; iat2 < sstr.atoms.size(); iat2++) {
      if (iat2 != iat1 && sstr.atoms[iat1].type == ispecie_min && sstr.atoms[iat2].type == ispecie_min) {
        fdisp = sstr.atoms[iat2].fpos - sstr.atoms[iat1].fpos;
        cdisp = sstr.atoms[iat2].cpos - sstr.atoms[iat1].cpos;
        if (IsTranslationFVector(sstr, fdisp)) {
          candidate_lattice_vector.push_back(cdisp);
        }
      }
    }
  }
  //  cerr << candidate_lattice_vector.size() << endl;
  int cnt = 0;
  // now generate triplets
  for (size_t iu = 0; iu < candidate_lattice_vector.size(); iu++) {
    for (uint i = 1; i <= 3; i++) {
      plattice(1, i) = candidate_lattice_vector.at(iu)(i);
    }
    for (size_t iv = 0; iv < candidate_lattice_vector.size() && iv != iu; iv++) {
      for (uint i = 1; i <= 3; i++) {
        plattice(2, i) = candidate_lattice_vector.at(iv)(i);
      }
      for (size_t iw = 0; iw < candidate_lattice_vector.size() && iw != iv && iw != iu; iw++) {
        for (uint i = 1; i <= 3; i++) {
          plattice(3, i) = candidate_lattice_vector.at(iw)(i);
        }
        if (det(plattice) > tolerance && det(plattice) < sstr_volume && det(plattice) < det(olattice)) { // well defined
          if (isdifferent(plattice, olattice, 0.0001)) {
            plattice = MinkowskiBasisReduction(plattice); // Minkowski first
            plattice = NiggliUnitCellForm(plattice); // Niggli Second
            if (isdifferent(plattice, olattice, 0.0001)) {
              olattice = plattice;
              //      cerr << det(olattice) << " " << cnt<< endl;
              cnt++;
            }
          }
        }
      }
    }
  }
  plattice = olattice;
  // done

  _atom atom;
  xstructure b = sstr;
  b.atoms.clear();
  b.lattice = plattice;
  b.lattice = roundoff(b.lattice, tolerance);
  b.FixLattices();
  b.write_lattice_flag = false;
  b.write_klattice_flag = false;
  b.write_DEBUG_flag = false;
  bool atom_found = false;
  // for scanning around
  //  double radius=aurostd::modulus(sstr.lattice(1))+aurostd::modulus(sstr.lattice(2))+aurostd::modulus(sstr.lattice(3));
  //  int dims=max(LatticeDimensionSphere(sstr.lattice,radius));
  const int dims = max(LatticeDimensionSphere(sstr.lattice, 1.5 * RadiusSphereLattice(sstr.lattice)));
  for (size_t iat = 0; iat < b.num_each_type.size(); iat++) {
    b.num_each_type[iat] = 0; // create enough space
    b.comp_each_type.at(iat) = 0; // create enough space
  }

  for (size_t iat = 0; iat < sstr.atoms.size(); iat++) {
    for (int i = -dims; i <= dims; i++) {
      for (int j = -dims; j <= dims; j++) {
        for (int k = -dims; k <= dims; k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
          atom = sstr.atoms[iat];
          atom.cpos = atom.cpos + ((double) i) * sstr.lattice(1) + ((double) j) * sstr.lattice(2) + ((double) k) * sstr.lattice(3);
          atom.fpos = C2F(b.lattice, atom.cpos);
          if (atom.fpos(1) >= -tolerance && atom.fpos(1) < 1.0 - tolerance && atom.fpos(2) >= -tolerance && atom.fpos(2) < 1.0 - tolerance && atom.fpos(3) >= -tolerance && atom.fpos(3) < 1.0 - tolerance) { // found something inside
            for (size_t ii = 0; ii < b.atoms.size() && !atom_found; ii++) {
              atom_found = identical(atom.cpos, b.atoms[ii].cpos, 0.1); // look in all the list of operations
            }
            // atom_found=false;
            if (!atom_found) {
              atom.fpos = roundoff(atom.fpos, tolerance);
              atom.cpos = roundoff(atom.cpos, tolerance);
              b.atoms.push_back(atom);
              b.num_each_type.at(atom.type)++; // CONVASP_MODE
              b.comp_each_type.at(atom.type) += atom.partial_occupation_value; // CONVASP_MODE
            }
          }
        }
      }
    }
  }
  b.GetStoich(); // CO20170724
  // rescale back to original scale.
  b.SetVolume(Volume(a) * b.atoms.size() / a.atoms.size());
  b = ReScale(b, a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b = LatticeReduction(b);
  // Put everything in new primitive cell.
  b = BringInCell(b);
  // check !
  const double fraction = Volume(a) / Volume(b);
  if (std::abs(b.atoms.size() * fraction - a.atoms.size()) > 0.1) {
    message << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    message << "        supercell has the wrong number of atoms" << endl;
    message << "        volume original    = " << Volume(a) << endl;
    message << "        volume prim        = " << Volume(b) << endl;
    message << "        a.scale            = " << a.scale << endl;
    message << "        b.scale            = " << b.scale << endl;
    message << "        a.atoms.size()     = " << a.atoms.size() << endl;
    message << "        b.atoms.size()     = " << b.atoms.size() << endl;
    message << "        fraction           = " << fraction << endl;
    message << "        supercell atoms    = " << fraction * b.atoms.size() << endl;
    message << b << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }
  // everything ok
  b.primitive_calculated = true; // DX20201007
  return b;
}

// second try
xstructure GetPrimitive2(const xstructure& a) {
  stringstream message;
  if (a.equiv_fpos_epsilon < 1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
  const double tolerance = a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);

  xstructure sstr = a;
  //  sstr.lattice=lattice;
  // sstr.scale=scale;
  std::vector<xvector<double>> candidate_lattice_vector;
  xmatrix<double> plattice(3, 3);
  vector<xmatrix<double>> plattice_list;
  const xvector<double> fdisp(3);
  xvector<double> cdisp(3);

  // Get all the data from the structure for easy use.
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr = ReScale(sstr, 1.0);
  // Put everything in the unit cell.
  //  cerr << sstr.scale << endl;
  sstr = BringInCell(sstr);
  // cerr << sstr.scale << endl;
  // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
  const string title = sstr.title;
  int i;
  for (size_t iat = 0; iat < sstr.num_each_type.size(); iat++) {
    sstr.num_each_type[iat] = 0; // create enough space
    sstr.comp_each_type.at(iat) = 0; // create enough space
  }
  // Create list of candidate plvs.

  candidate_lattice_vector.push_back(sstr.lattice(1));
  candidate_lattice_vector.push_back(sstr.lattice(2));
  candidate_lattice_vector.push_back(sstr.lattice(3));

  for (size_t ia = 1; ia < sstr.atoms.size(); ia++) {
    // Here we calculate all displacements from the first atom, i.e., the
    // first atom is considered the origin.
    cdisp = sstr.atoms[ia].cpos - sstr.atoms.at(0).cpos;
    if (IsTranslationCVector(sstr, cdisp)) {
      candidate_lattice_vector.push_back(cdisp);
    }
  }
  // Now the candidate_lattice_vector have been found and we must take all possible
  // traids of *distinct* vectors and get the smallest volume.
  const int num_cand = candidate_lattice_vector.size();
  double min_vol;
  double new_vol;
  min_vol = GetVol(sstr.lattice);
  for (int iu = 0; iu < num_cand; iu++) {
    for (i = 1; i <= 3; i++) {
      plattice(1, i) = candidate_lattice_vector.at(iu)(i);
    }
    for (int iv = 0; iv < num_cand; iv++) {
      for (i = 1; i <= 3; i++) {
        plattice(2, i) = candidate_lattice_vector.at(iv)(i);
      }
      for (int iw = 0; iw < num_cand; iw++) {
        if (iu != iw && iu != iv && iv != iw) {
          for (i = 1; i <= 3; i++) {
            plattice(3, i) = candidate_lattice_vector.at(iw)(i);
          }
          new_vol = GetVol(plattice);
          // if(new_vol<=min_vol && new_vol>tolerance) min_vol=new_vol;
          if (new_vol <= min_vol && new_vol > 1.0e-5) {
            min_vol = new_vol;
          }
        }
      } // iw loop
    } // iv loop
  } // iu loop
  // Now that we know the min_vol we must go back through all the triads
  // and store the set that have the min_vol.
  for (int iu = 0; iu < num_cand; iu++) {
    for (i = 1; i <= 3; i++) {
      plattice(1, i) = candidate_lattice_vector.at(iu)(i);
    }
    for (int iv = 0; iv < num_cand; iv++) {
      for (i = 1; i <= 3; i++) {
        plattice(2, i) = candidate_lattice_vector.at(iv)(i);
      }
      for (int iw = 0; iw < num_cand; iw++) {
        if (iu != iw && iu != iv && iv != iw) {
          for (i = 1; i <= 3; i++) {
            plattice(3, i) = candidate_lattice_vector.at(iw)(i);
          }
          new_vol = GetVol(plattice);
          if (aurostd::abs(new_vol - min_vol) < tolerance && new_vol > 1.0e-5) {
            plattice_list.push_back(plattice); // Add to plattice_list
          }
          // if new_vol=min_vol
        } // if iu!=iv!=iw
      } // iw loop
    } // iv loop
  } // iu loop
  // Now we have a set of all possible primitive lattices
  // and we want to pick one that looks as much like the
  // original lattice as possible.
  double paw;
  double pawopt = 1e6;
  for (size_t i = 0; i < plattice_list.size(); i++) {
    paw = 0;
    for (int j = 1; j <= 3; j++) {
      paw += angle(plattice_list.at(i)(j), sstr.lattice(j));
    }
    if (paw < pawopt) {
      pawopt = paw;
      plattice = plattice_list[i];
    }; // minimization
  }
  // For now just pick the first one.
  // if(plattice_list.size()>0) plattice=plattice_list.at(1);
  // If volume did not reduce then just keep original lattice
  if (aurostd::abs(GetVol(sstr.lattice) - GetVol(plattice)) < tolerance) {
    plattice = sstr.lattice;
  }

  // Now we have the plvs for the new structure.  Now we must
  // construct a new structure object with the plvs and return it.
  // This new structure object will be just like the old one but
  // with new lattice vectors and new basis atoms (which will be a
  // subset of the original ones, chosen to have direct coordinates
  // less than 1).

  // To create new structure we must set all the parameters.
  // To get new basis atoms, loop over all the original basis atoms,
  // get new direct coordinates, make all less than one, then take unique
  // values.

  _atom atom;
  xstructure b = sstr;
  b.atoms.clear();
  b.lattice = plattice;
  b.lattice = roundoff(b.lattice, tolerance);
  b.FixLattices();
  b.write_lattice_flag = true;
  b.write_klattice_flag = false;
  b.write_DEBUG_flag = false;
  bool atom_found = false;
  // for scanning around
  double radius;
  radius = aurostd::modulus(sstr.lattice(1)) + aurostd::modulus(sstr.lattice(2)) + aurostd::modulus(sstr.lattice(3));
  const int dims = max(LatticeDimensionSphere(sstr.lattice, radius));
  for (size_t iat = 0; iat < b.num_each_type.size(); iat++) {
    b.num_each_type[iat] = 0; // create enough space
    b.comp_each_type.at(iat) = 0; // create enough space
  }

  for (size_t iat = 0; iat < sstr.atoms.size(); iat++) {
    for (int i = -dims; i <= dims; i++) {
      for (int j = -dims; j <= dims; j++) {
        for (int k = -dims; k <= dims; k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
          atom = sstr.atoms[iat];
          atom.cpos = atom.cpos + ((double) i) * sstr.lattice(1) + ((double) j) * sstr.lattice(2) + ((double) k) * sstr.lattice(3);
          atom.fpos = C2F(b.lattice, atom.cpos);
          if (atom.fpos(1) >= -tolerance && atom.fpos(1) < 1.0 - tolerance && atom.fpos(2) >= -tolerance && atom.fpos(2) < 1.0 - tolerance && atom.fpos(3) >= -tolerance && atom.fpos(3) < 1.0 - tolerance) { // found something inside
            for (size_t ii = 0; ii < b.atoms.size() && !atom_found; ii++) {
              atom_found = identical(atom.cpos, b.atoms[ii].cpos, 0.1); // look in all the list of operations
            }
            // atom_found=false;
            if (!atom_found) {
              atom.fpos = roundoff(atom.fpos, tolerance);
              atom.cpos = roundoff(atom.cpos, tolerance);
              b.atoms.push_back(atom);
              b.num_each_type.at(atom.type)++; // CONVASP_MODE
              b.comp_each_type.at(atom.type) += atom.partial_occupation_value; // CONVASP_MODE
            }
          }
        }
      }
    }
  }
  b.GetStoich(); // CO20170724
  // some check
  const double fraction = (GetVol(sstr.lattice) / GetVol(b.lattice));
  if (std::abs(b.atoms.size() * fraction - a.atoms.size()) > 0.1) {
    cerr << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    cerr << "        supercell has the wrong number of atoms" << endl;
    cerr << "        volume original    = " << GetVol(sstr.lattice) << endl;
    cerr << "        volume prim        = " << GetVol(b.lattice) << endl;
    // cerr << "        sstr.scale            = " << sstr.scale << endl;
    // cerr << "        b.scale            = " << b.scale << endl;
    cerr << "        sstr.atoms.size()     = " << a.atoms.size() << endl;
    cerr << "        b.atoms.size()     = " << b.atoms.size() << endl;
    cerr << "        fraction           = " << fraction << endl;
    cerr << "        supercell atoms    = " << fraction * b.atoms.size() << endl;
    cerr << GetVol(b.lattice) << endl;
    cerr << GetVol(a.lattice) << endl;
    //   throw aurostd::xerror(__AFLOW_FILE__,XPID+"xstructure::GetPrimitive():","Throw for debugging purposes.",_GENERIC_ERROR_);
  }
  // rescale back to original scale.
  b = ReScale(b, a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b = LatticeReduction(b);
  // Put everything in new primitive cell.
  b = BringInCell(b);
  // everything ok
  b.primitive_calculated = true; // DX20201007
  return b;
}

// third try
xstructure GetPrimitive3(const xstructure& a) {
  stringstream message;
  if (a.equiv_fpos_epsilon < 1.0e-12) {
    message << "Zero tolerance: " << a.equiv_fpos_epsilon;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }
  const double tolerance = a.equiv_fpos_epsilon;

  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);
  xstructure sstr = a;
  sstr.SetVolume(sstr.atoms.size());
  sstr = ReScale(sstr, 1.0);
  sstr = BringInCell(sstr);

  const double sstr_volume = sstr.Volume();
  _aflags aflags;
  aflags.Directory = "./";
  // identify the minimum set of atoms
  const bool PGROUPWRITE = true;
  const bool FGROUPWRITE = true; //,IATOMSWRITE=true;
  const bool OSSWRITE = false; // to FileMESSAGE, does not matter as it is /dev/null
  ofstream FileMESSAGE("/dev/null");
  const bool _QUIET_ = XHOST.QUIET;
  XHOST.QUIET = true;
  SYM::CalculatePointGroup(FileMESSAGE, sstr, aflags, PGROUPWRITE, OSSWRITE, cout);
  SYM::CalculateFactorGroup(FileMESSAGE, sstr, aflags, FGROUPWRITE, OSSWRITE, cout);
  //  SYM::CalculateInequivalentAtoms(FileMESSAGE,sstr,aflags,IATOMSWRITE,OSSWRITE,cout);

  xmatrix<double> plattice(3, 3);
  xmatrix<double> olattice(3, 3);
  const xvector<double> fdisp(3);
  const xvector<double> cdisp(3);
  std::vector<xvector<double>> candidate_lattice_vector;

  // generate list of vectors
  candidate_lattice_vector.push_back(sstr.lattice(1)); // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(2)); // lattice is made of good vectors
  candidate_lattice_vector.push_back(sstr.lattice(3)); // lattice is made of good vectors
  olattice = sstr.lattice; // the lattice is always a good lattice

  bool sym_found;
  for (size_t i = 0; i < sstr.fgroup.size(); i++) {
    if (aurostd::modulus(sstr.fgroup[i].ctau) > 0.01) {
      //    cerr << i << " " << sstr.fgroup[i].ctau << endl;
      sym_found = false;
      for (size_t ii = 0; ii < candidate_lattice_vector.size() && !sym_found; ii++) {
        sym_found = identical(sstr.fgroup[i].ctau, candidate_lattice_vector[ii], tolerance);
      }
      // look in all the list of operations
      if (sym_found == false) { // new operation, generate and save it
        candidate_lattice_vector.push_back(sstr.fgroup[i].ctau);
        cerr << i << " " << sstr.fgroup[i].ctau << endl;
      }
    }
  }

  cerr << candidate_lattice_vector.size() << endl;
  int cnt = 0;
  // now generate triplets
  for (size_t iu = 0; iu < candidate_lattice_vector.size(); iu++) {
    for (uint i = 1; i <= 3; i++) {
      plattice(1, i) = candidate_lattice_vector.at(iu)(i);
    }
    for (size_t iv = 0; iv < candidate_lattice_vector.size() && iv != iu; iv++) {
      for (uint i = 1; i <= 3; i++) {
        plattice(2, i) = candidate_lattice_vector.at(iv)(i);
      }
      for (size_t iw = 0; iw < candidate_lattice_vector.size() && iw != iv && iw != iu; iw++) {
        for (uint i = 1; i <= 3; i++) {
          plattice(3, i) = candidate_lattice_vector.at(iw)(i);
        }
        if (det(plattice) > tolerance && det(plattice) < sstr_volume && det(plattice) < det(olattice)) { // well defined
          if (isdifferent(plattice, olattice, 0.0001)) {
            plattice = MinkowskiBasisReduction(plattice); // Minkowski first
            plattice = NiggliUnitCellForm(plattice); // Niggli Second
            if (isdifferent(plattice, olattice, 0.0001)) {
              olattice = plattice;
              //      cerr << det(olattice) << " " << cnt<< endl;
              cnt++;
            }
          }
        }
      }
    }
  }
  plattice = olattice;
  // done

  _atom atom;
  xstructure b = sstr;
  b.atoms.clear();
  b.lattice = plattice;
  b.lattice = roundoff(b.lattice, tolerance);
  b.FixLattices();
  b.write_lattice_flag = false;
  b.write_klattice_flag = false;
  b.write_DEBUG_flag = false;
  bool atom_found = false;
  // for scanning around
  const int dims = max(LatticeDimensionSphere(sstr.lattice, RadiusSphereLattice(sstr.lattice)));
  for (size_t iat = 0; iat < b.num_each_type.size(); iat++) {
    b.num_each_type[iat] = 0; // create enough space
    b.comp_each_type.at(iat) = 0; // create enough space
  }

  for (size_t iat = 0; iat < sstr.atoms.size(); iat++) {
    for (int i = -dims; i <= dims; i++) {
      for (int j = -dims; j <= dims; j++) {
        for (int k = -dims; k <= dims; k++) {
          //    atom=BringInCell(sstr.atoms.at(i),sstr.lattice);
          atom = sstr.atoms[iat];
          atom.cpos = atom.cpos + ((double) i) * sstr.lattice(1) + ((double) j) * sstr.lattice(2) + ((double) k) * sstr.lattice(3);
          atom.fpos = C2F(b.lattice, atom.cpos);
          if (atom.fpos(1) >= -tolerance && atom.fpos(1) < 1.0 - tolerance && atom.fpos(2) >= -tolerance && atom.fpos(2) < 1.0 - tolerance && atom.fpos(3) >= -tolerance && atom.fpos(3) < 1.0 - tolerance) { // found something inside
            for (size_t ii = 0; ii < b.atoms.size() && !atom_found; ii++) {
              atom_found = identical(atom.cpos, b.atoms[ii].cpos, 0.1); // look in all the list of operations
            }
            // atom_found=false;
            if (!atom_found) {
              atom.fpos = roundoff(atom.fpos, tolerance);
              atom.cpos = roundoff(atom.cpos, tolerance);
              b.atoms.push_back(atom);
              b.num_each_type.at(atom.type)++; // CONVASP_MODE
              b.comp_each_type.at(atom.type) += atom.partial_occupation_value; // CONVASP_MODE
            }
          }
        }
      }
    }
  }
  b.GetStoich(); // CO20170724
  // rescale back to original scale.
  b.SetVolume(Volume(a) * b.atoms.size() / a.atoms.size());
  b = ReScale(b, a.scale);
  // fix it up with the new Minkowsky and Niggli reductions
  b = LatticeReduction(b);
  XHOST.QUIET = _QUIET_;

  // Put everything in new primitive cell.
  b = BringInCell(b);
  // check !
  const double fraction = Volume(a) / Volume(b);
  if (std::abs(b.atoms.size() * fraction - a.atoms.size()) > 0.1) {
    message << "ERROR   xstructure xstructure::GetPrimitive(void)" << endl;
    message << "        supercell has the wrong number of atoms" << endl;
    message << "        volume original    = " << Volume(a) << endl;
    message << "        volume prim        = " << Volume(b) << endl;
    message << "        a.scale            = " << a.scale << endl;
    message << "        b.scale            = " << b.scale << endl;
    message << "        a.atoms.size()     = " << a.atoms.size() << endl;
    message << "        b.atoms.size()     = " << b.atoms.size() << endl;
    message << "        fraction           = " << fraction << endl;
    message << "        supercell atoms    = " << fraction * b.atoms.size() << endl;
    message << b << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }
  // everything ok
  b.primitive_calculated = true; // DX20201007
  return b;
}

// ***************************************************************************
// Operator GetPrimitive
// ***************************************************************************

void xstructure::GetPrimitive() {
  extern xstructure GetPrimitive(const xstructure& a); // so it does not recurse
  const xstructure a(*this);
  *this = GetPrimitive(a);
}

void xstructure::GetPrimitive(double tolerance) {
  extern xstructure GetPrimitive(const xstructure& a, double tolerance); // so it does not recurse
  const xstructure a(*this);
  *this = GetPrimitive(a, tolerance);
}

void xstructure::GetPrimitive2() {
  extern xstructure GetPrimitive2(const xstructure& a); // so it does not recurse
  const xstructure a(*this);
  *this = GetPrimitive2(a);
}

void xstructure::GetPrimitive3() {
  extern xstructure GetPrimitive3(const xstructure& a); // so it does not recurse
  const xstructure a(*this);
  *this = GetPrimitive3(a);
}

// ***************************************************************************
// Function GetPrimitiveStructures() //DX20201006
// ***************************************************************************
void GetPrimitiveStructures(vector<xstructure>& structures, uint start_index, uint end_index) {
  // Converts a set of xstructures to their primitive representation
  // Optional indices can be included; useful for pre-distributed
  // threading schemes
  // Default: run over entire range

  // if end index is greater than structures.size(), then compute primitive cell for all structures
  if (end_index > structures.size()) {
    end_index = structures.size();
  }

  for (uint i = start_index; i < end_index; i++) {
    structures[i].GetPrimitive();
  }
}

// ***************************************************************************
// Function MinDist
// ***************************************************************************
double xstructure::MinDist() {
  dist_nn_min = SYM::minimumDistance(*this);
  return dist_nn_min;
}

// ***************************************************************************
// Function NearestNeighbor() // moved from aflow_xproto.cpp
// ***************************************************************************
double NearestNeighbor(const xstructure& str_in) {
  return SYM::minimumDistance(str_in);
}

// ***************************************************************************
// Function NearestNeighbors() //DX20201230 - moved from XtalFinder
// ***************************************************************************
vector<double> NearestNeighbors(const xstructure& xstr) {
  // Determine the nearest neighbor distances centered on each atom
  // of the structure (needed for XtalFinder)

  vector<double> all_nn_distances;
  double nn = AUROSTD_MAX_DOUBLE;

  for (size_t i = 0; i < xstr.atoms.size(); i++) {
    nn = NearestNeighborToAtom(xstr, i);
    all_nn_distances.push_back(nn);
  }
  return all_nn_distances;
}

// ***************************************************************************
// Function NearestNeighborToAtom() //DX20201230 - moved from XtalFinder
// ***************************************************************************
double NearestNeighborToAtom(const xstructure& xstr, uint k) {
  // Find the minimum interatomic distance in the structure to atom k
  // Different than SYM::minimumDistance(): only considers one atom index
  // in minimization routine, as opposed to global minimimum
  // (considering one atom only affords speed ups)
  // Use resetLatticeDimension() to update search radius for nearest
  // neighbors: once we find a neighbor, update/reduce how far we
  // need to search to find a closer neighbor

  double min_dist = AUROSTD_MAX_DOUBLE;
  double prev_min_dist = 0; // DX20190716
  const xmatrix<double> lattice = xstr.lattice; // NEW

  // DX speed increase
  // perhaps can speed up even more, since the lattice doesn't change for the xstr...
  vector<xvector<double>> l1;
  vector<xvector<double>> l2;
  vector<xvector<double>> l3;
  vector<int> a_index;
  vector<int> b_index;
  vector<int> c_index;
  xvector<int> dims(3); // DX20190710 - use robust method
  dims[1] = dims[2] = dims[3] = 0; // reset

  xvector<double> tmp_coord;
  xvector<double> incell_dist;
  xvector<double> a_component;
  xvector<double> ab_component; // DX20200329
  double incell_mod = AUROSTD_MAX_DOUBLE;

  uint ii = 0;
  uint m = 0;
  uint n = 0;
  uint p = 0;
  uint m_size = 0;
  uint n_size = 0;
  uint p_size = 0;

  for (ii = 0; ii < xstr.atoms.size(); ii++) {
    if (ii != k) {
      if (min_dist < prev_min_dist) {
        if (!(dims[1] == 1 && dims[2] == 1 && dims[3] == 1)) {
          // update the dimensions based on new search radius (min_dist)
          resetLatticeDimensions(lattice, min_dist, dims, l1, l2, l3, a_index, b_index, c_index);
          prev_min_dist = min_dist;
          m_size = l1.size();
          n_size = l2.size();
          p_size = l3.size();
        }
      }
      incell_dist = xstr.atoms[k].cpos - xstr.atoms[ii].cpos;
      incell_mod = aurostd::modulus(incell_dist);
      if (incell_mod < min_dist) {
        if (!(dims[1] == 1 && dims[2] == 1 && dims[3] == 1)) {
          // update the dimensions based on new search radius (incell_mod)
          resetLatticeDimensions(lattice, incell_mod, dims, l1, l2, l3, a_index, b_index, c_index);
          m_size = l1.size();
          n_size = l2.size();
          p_size = l3.size();
        }
        prev_min_dist = incell_mod;
      }
      // DX20180423 - running vector in each loop saves computations; fewer duplicate operations
      for (m = 0; m < m_size; m++) {
        a_component = incell_dist + l1[m]; // DX : coord1-coord2+a*lattice(1)
        for (n = 0; n < n_size; n++) {
          ab_component = a_component + l2[n]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
          for (p = 0; p < p_size; p++) {
            tmp_coord = ab_component + l3[p]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
            min_dist = aurostd::min(min_dist, aurostd::modulus(tmp_coord));
          }
        }
      }
    }
  }

  return min_dist;
}

// ***************************************************************************
// Function ReScale
// ***************************************************************************
xstructure ReScale(const xstructure& a, const double& in_scale) {
  // This resets scale and changes the cell parameters and coordinates
  // appropriately.  Keeps volume fixed.
  if (in_scale == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "ReScale()", "in_scale must be non zero", _INPUT_ILLEGAL_);
  }
  xstructure b(a);
  if (aurostd::identical(b.scale, in_scale, _ZERO_TOL_)) {
    return b;
  }
  // try hard not to introduce precision errors, currently we print scale with precision 6
  b.lattice = b.lattice * b.scale / in_scale;
  b.origin = b.origin * b.scale / in_scale;
  b.f2c = trasp(b.lattice);
  b.c2f = inverse(trasp(b.lattice));
  // klattice already contained the scale so it does not need to be fixed
  // b.klattice=b.klattice*in_scale/b.scale;
  for (int i = 0; i < (int) b.atoms.size(); i++) {
    b.atoms[i].fpos = a.atoms.at(i).fpos;
    b.atoms[i].cpos = a.atoms.at(i).cpos * b.scale / in_scale;
  }
  if (b.fgroup_calculated) {
    for (int fg = 0; fg < (int) b.fgroup.size(); fg++) {
      b.fgroup[fg].ctau = b.fgroup[fg].ctau * b.scale / in_scale;
    }
  }
  if (b.sgroup_calculated) {
    for (int sg = 0; sg < (int) b.sgroup.size(); sg++) {
      b.sgroup[sg].ctau = b.sgroup[sg].ctau * b.scale / in_scale;
      b.sgroup[sg].ctrasl = b.sgroup[sg].ctrasl * b.scale / in_scale;
    }
  }
  b.scale = in_scale;
  b.FixLattices(); // touched scale, then fix the lattices
  return b;
}

void xstructure::ReScale(const double& in_scale) {
  if (in_scale == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "ReScale()", "in_scale must be non zero", _INPUT_ILLEGAL_);
  }
  if (aurostd::identical(scale, in_scale, _ZERO_TOL_)) {
    return;
  }
  // try hard not to introduce precision errors, currently we print scale with precision 6
  lattice = lattice * scale / in_scale;
  origin = origin * scale / in_scale;
  f2c = trasp(lattice);
  c2f = inverse(trasp(lattice));
  // klattice already contained the scale so it does not need to be fixed
  // klattice=klattice*in_scale/scale;
  for (int i = 0; i < (int) atoms.size(); i++) {
    atoms[i].cpos = atoms[i].cpos * scale / in_scale;
  }
  if (fgroup_calculated) {
    for (int fg = 0; fg < (int) fgroup.size(); fg++) {
      fgroup[fg].ctau = fgroup[fg].ctau * scale / in_scale;
    }
  }
  if (sgroup_calculated) {
    for (int sg = 0; sg < (int) sgroup.size(); sg++) {
      sgroup[sg].ctau = sgroup[sg].ctau * scale / in_scale;
      sgroup[sg].ctrasl = sgroup[sg].ctrasl * scale / in_scale;
    }
  }
  scale = in_scale;
  FixLattices(); // touched scale, then fix the lattices
}

// ***************************************************************************
// Function SetScale
// ***************************************************************************
void xstructure::SetScale(const double& in_scale) {
  if (in_scale == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "SetScale()", "in_scale must be non zero", _INPUT_ILLEGAL_);
  }
  scale = in_scale;
  FixLattices(); // touched scale, then fix the lattices
}

xstructure SetScale(const xstructure& a, const double& in_scale) {
  // This resets scale.  Keeps volume fixed.
  if (in_scale == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "SetScale()", "in_scale must be non zero", _INPUT_ILLEGAL_);
  }
  xstructure b;
  b = a;
  b.scale = in_scale;
  b.FixLattices(); // touched scale, then fix the lattices
  return b;
}

// AS20200514 START
//  ***************************************************************************
//  Function UpdateCartesianCoordinates
//  ***************************************************************************
void xstructure::UpdateCartesianCoordinates() {
  for (size_t at = 0; at < atoms.size(); at++) {
    atoms[at].cpos = scale * f2c * atoms[at].fpos;
  }
}

// AS20200514 END

// ***************************************************************************
// Function SetVolume
// ***************************************************************************
void xstructure::SetVolume(const double& in_volume) {
  if (in_volume == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, "SetVolume()", "in_scale must be non zero", _INPUT_ILLEGAL_);
  } // CO20200201
  if (det(lattice) < 0.0) { // CO20200201
    stringstream message; // CO20200201
    message << "Found negative determinant for lattice (det()=" << det(lattice) << "). Flip your basis."; // CO20200201
    throw aurostd::xerror(__AFLOW_FILE__, "SetVolume()", message, _INPUT_ILLEGAL_); // CO20200201
  } // CO20200201
  scale = std::pow((double) in_volume / det(lattice), 1.0 / 3.0);
  FixLattices(); // touched scale, then fix the lattices
}

xstructure SetVolume(const xstructure& a, const double& in_volume) {
  xstructure b(a);
  b.SetVolume(in_volume);
  return b;
}

// ***************************************************************************
// Function SetAutoVolume
// ***************************************************************************
void xstructure::SetAutoVolume(bool use_AFLOW_defaults_in) { // CO20191010
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " xstr=" << endl << (*this) << endl;
    cerr << __AFLOW_FUNC__ << " xstr.species=" << aurostd::joinWDelimiter((*this).species, ",") << endl;
    cerr << __AFLOW_FUNC__ << " fixing volume" << endl;
    cerr << __AFLOW_FUNC__ << " volume_orig=" << GetVolume() << endl;
  }
  double volume = 0; //,voli=0;
  bool use_AFLOW_defaults = use_AFLOW_defaults_in;
  // try and pull from species_volume first
  for (size_t i = 0; i < atoms.size() && !use_AFLOW_defaults; i++) {
    for (size_t j = 0; j < num_each_type.size() && !use_AFLOW_defaults; j++) {
      if (atoms[i].name == species[j]) {
        const double& voli = species_volume[j];
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " atoms[i].name=" << atoms[i].name << " atoms[i].vol=" << voli << endl;
        }
        if (voli == NNN || aurostd::isequal(voli, 0.0, _ZERO_TOL_)) {
          use_AFLOW_defaults = true;
        }
        if (aurostd::isequal(atoms[i].partial_occupation_value, 0.0, _ZERO_TOL_)) {
          message << "partial_occupation_value==0.0 for atom=" << i;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        }
        volume += atoms[i].partial_occupation_value * voli;
      }
    }
  }
  // otherwise get defaults from AFLOW
  if (use_AFLOW_defaults || std::abs(volume) < _XPROTO_ZERO_VOL_) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " using automatic volumes" << endl;
    }
    volume = 0;
    double voli = 0.0;
    for (size_t i = 0; i < atoms.size(); i++) {
      for (size_t j = 0; j < num_each_type.size(); j++) {
        if (atoms[i].name == species[j]) {
          voli = GetAtomVolume(atoms[i].name);
          if (voli == NNN || aurostd::isequal(voli, 0.0, _ZERO_TOL_)) {
            message << "No volume found for " << atoms[i].name << " (auto volumes)";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
          }
          if (aurostd::isequal(atoms[i].partial_occupation_value, 0.0, _ZERO_TOL_)) {
            message << "partial_occupation_value==0.0 for atom=" << i;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
          }
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " atoms[i=" << i << "].name=" << atoms[i].name << " pocc=" << atoms[i].partial_occupation_value << " vol=" << voli << endl;
          }
          volume += atoms[i].partial_occupation_value * voli;
        }
      }
    }
  }
  if (std::abs(volume) < _XPROTO_ZERO_VOL_) {
    message << "Final volume==0, check species default volumes";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
  }
  SetVolume(volume);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " volume_new=" << GetVolume() << endl;
  }
}

// ***************************************************************************
// Function InflateLattice
// ***************************************************************************
void xstructure::InflateLattice(const double& coefficient) {
  if (coefficient == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::InflateLattice()", "coefficient must be non zero.", _INPUT_ILLEGAL_);
  }
  //  scale=coefficient*scale;
  lattice = coefficient * lattice;
  FixLattices(); // touched scale/lattice, then fix the lattices
}

xstructure InflateLattice(const xstructure& a, const double& coefficient) {
  // This resets scale.  Keeps volume fixed.
  if (coefficient == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::InflateLattice()", "coefficient must be non zero.", _INPUT_ILLEGAL_);
  }
  xstructure b;
  b = a;
  //  b.scale=coefficient*b.scale;
  b.lattice = coefficient * b.lattice;
  b.FixLattices(); // touched scale/lattice, then fix the lattices
  return b;
}

// ***************************************************************************
// Function InflateVolume
// ***************************************************************************
void xstructure::InflateVolume(const double& coefficient) {
  if (coefficient == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::InflateVolume()", "coefficient must be non zero", _INPUT_ILLEGAL_);
  }
  // scale=std::pow((double) coefficient,(double) 1/3)*scale;
  lattice = std::pow((double) coefficient, (double) 1 / 3) * lattice;
  FixLattices(); // touched scale/lattice, then fix the lattices
  UpdateCartesianCoordinates(); // AS20200514
}

xstructure InflateVolume(const xstructure& a, const double& coefficient) {
  if (coefficient == 0.0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::InflateVolume()", "coefficient must be non zero", _INPUT_ILLEGAL_);
  }
  xstructure b;
  b = a;
  //  b.scale=std::pow((double) coefficient,(double) 1/3)*b.scale;
  b.lattice = std::pow((double) coefficient, (double) 1 / 3) * b.lattice;
  b.FixLattices(); // touched scale/lattice, need to fix the lattices
  b.UpdateCartesianCoordinates(); // AS20200514
  return b;
}

// ***************************************************************************
// Function GetVolume
// ***************************************************************************
double xstructure::GetVolume() const { // CO20200201
  return scale * scale * scale * det(lattice);
}

double GetVolume(const xstructure& a) {
  return a.GetVolume();
}

// ***************************************************************************
// Function Volume
// ***************************************************************************
double xstructure::Volume() const { // CO20200201
  return GetVolume(); // CO20200201
  //[CO20200201]return scale*scale*scale*det(lattice);
}

double Volume(const xstructure& a) {
  return a.GetVolume(); // CO20200201
  //[CO20200201]return a.scale*a.scale*a.scale*det(a.lattice);
}

_atom BringCloseToOrigin(_atom& atom, xmatrix<double>& f2c) {
  _atom atom_out = atom;
  const xvector<double> v_in = atom.fpos;
  const xvector<int> ijk = atom.ijk;
  for (uint i = 1; i < 4; i++) {
    while (v_in(i) > (1.0 - _ZERO_TOL_)) {
      v_in(i) = v_in(i) - 1;
      ijk(i) -= 1;
    }
    while (v_in(i) <= -_ZERO_TOL_) { // fixed DX
      v_in(i) = v_in(i) + 1;
      ijk(i) += 1;
    }
  }
  atom_out.fpos = v_in;
  atom_out.cpos = f2c * atom_out.fpos;
  atom_out.ijk = ijk;
  return atom_out;
}

bool uniqueAtomInCell(_atom& atom, deque<_atom>& atoms) {
  if (inCell(atom.fpos)) {
    if (!alreadyInCell(atom, atoms)) {
      return true;
    }
    return false;
  }
  return false;
}

// DX START
// DX END

// ***************************************************************************
// atomInCell()
// ***************************************************************************
bool atomInCell(const _atom& atom, double tolerance, double upper_bound, double lower_bound) {
  // ME+DX20210203 - added bounds

  // check if the atom is in the unit cell based on fractional coordinates
  // if you use the non-default tolerance (i.e., _ZERO_TOL_), this alone is not robust
  // Note: check over each component and returning false immediately (faster)

  return inCell(atom.fpos, tolerance, upper_bound, lower_bound);
}

// ***************************************************************************
// inCell()
// ***************************************************************************
// ME20210128 - Added bounds
bool inCell(const xvector<double>& pos_vec, double tolerance, double upper_bound, double lower_bound) {
  // check if the position is in the unit cell based on fractional coordinates
  // if you use the non-default tolerance (i.e., _ZERO_TOL_), this alone is not robust
  // Note: check over each component and returning false immediately (faster)

  for (uint f = 1; f < 4; f++) {
    // ME20210128: Used to be pos_vec[f] > 1.0 + tolerance.
    // Adjusted to use the same cut-off criterion as bringInCell
    if ((pos_vec[f] - upper_bound) >= -tolerance || (pos_vec[f] - lower_bound) < -tolerance) { // allows tunable cutoff
      return false;
    }
  }
  return true;
}

// DX20180726 - check if already in cell - START
bool alreadyInCell(_atom& atom, deque<_atom> atoms) {
  for (size_t i = 0; i < atoms.size(); i++) {
    xvector<double> fdiff = atom.fpos - atoms[i].fpos;
    fdiff = SYM::minimizeDistanceFractionalMethod(fdiff); // DX20190613
    if (aurostd::abs(fdiff(1)) < _ZERO_TOL_ && aurostd::abs(fdiff(2)) < _ZERO_TOL_ && aurostd::abs(fdiff(3)) < _ZERO_TOL_) {
      return true;
    }
  }
  return false;
}

// DX20180726 - check if already in cell - END

// ***************************************************************************
// Function GetSuperCell
// ***************************************************************************
// This funtion creates SuperCells with 9 or 3 elements...
// the old routine by Dane Morgan does not work, so I rewrote this
// one from scratch.
// The algorithm is simple: make the bigger cell
// b.lattice=supercell*a.lattice;
// and generate a bunch of atoms around the old cell and check
// if they are inside the new cell... There is a checksum with ERROR
// if the numbers do not match. Stefano Curtarolo (aug07).
// CO (2017) adding maps to/from supercell/primitive structure, get_symmetry and get_full_basis flags
// get_symmetry will propagate symmetry of primitive cell
// WARNING: if the structure is a derivative structure (i.e., non-uniform expansion, see derivative_structure) then
// NOT all symmetry operations will work, as the symmetry is reduced
// YOU need to check which of these don't work, unless you also calculate get_full_basis,
// which checks the validity of the symmetry operation by finding the full basis map

// CO START
xstructure GetSuperCell(const xstructure& aa, const xmatrix<double>& supercell, vector<int>& sc2pcMap, vector<int>& pc2scMap, bool get_symmetry, bool get_full_basis, bool force_supercell_matrix, bool force_strict_pc2scMap)
// DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
// xstructure GetSuperCell(const xstructure& aa, const xmatrix<double> &supercell)
{ // CO20200106 - patching for auto-indenting
  // #define _eps_scell_ 1.0e-5
  // #define _eps_scell_ 0.0
  //  check for error
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;
  const double vol_supercell = det(supercell);
  if (std::abs(vol_supercell) < 0.001) {
    message << "Singular supercell matrix";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }
  xstructure a(aa);
  a.ReScale(1.0); // the nuclear option, the only way not to mess around with scale EVERYWHERE
  // DO NOT MODIFY STRUCTURE IN HERE, WE WANT TO PROPAGATE SYMMETRY FROM PRIMITIVE STRUCTURE!
  // a.BringInCell();
  xstructure b(a); // actual supercell, need to copy!
  // pflow::CalculateFullSymmetry(a);
  _atom atom;
  const _atom atom2;
  _sym_op pSymOp;
  _sym_op fSymOp;
  _sym_op aSymOp; // CO
  int i;
  int j;
  int k; //,dim;
  uint iat; //,at_indx; //CO
  xvector<int> dims(3);
  xvector<double> cshift(3); // CO
  vector<xvector<double>> cshifts; // CO
  // double zeroTol=1e-10;
  double radius;
  b.lattice = supercell * a.lattice;
  // the scale is kept the same....it is always saved as positive
  b.FixLattices();
  for (i = 0; i < (int) b.num_each_type.size(); i++) {
    b.num_each_type[i] = 0;
  }
  for (i = 0; i < (int) b.comp_each_type.size(); i++) {
    b.comp_each_type[i] = 0;
  }
  b.atoms.clear();
  sc2pcMap.clear();
  pc2scMap.clear(); // CO20170722 - clear this everytime

  bool skew = false; // DX20190319 - declared outside loop

  if (get_symmetry) { // DX20190319 - added if statement; don't calc unless necessary
    if (b.dist_nn_min == AUROSTD_NAN) {
      b.dist_nn_min = SYM::minimumDistance(a);
    } // calculate dist_nn_min
    if (b.sym_eps == AUROSTD_NAN) {
      b.sym_eps = SYM::defaultTolerance(b);
    }
    skew = SYM::isLatticeSkewed(b.lattice, b.dist_nn_min, b.sym_eps); // DX20190319 - declared above
  } // DX //DX20190319 - added if statement; don't calc unless necessary

  const double nx = supercell(1, 1);
  const double ny = supercell(2, 2);
  const double nz = supercell(3, 3);
  // can only propagate if nx==ny==nz and diagonal
  // false means we have a derivative structure
  // bool derivative_structure=!(abs(supercell(1,1)-supercell(2,2))<zeroTol && abs(supercell(1,1)-supercell(3,3))<zeroTol &&
  //                             abs(supercell(1,2))<zeroTol && abs(supercell(1,3))<zeroTol &&
  //                             abs(supercell(2,1))<zeroTol && abs(supercell(2,3))<zeroTol &&
  //                             abs(supercell(3,1))<zeroTol && abs(supercell(3,2))<zeroTol);
  const bool derivative_structure = !(aurostd::isdiagonal(supercell) && aurostd::isequal(nx, ny) && aurostd::isequal(ny, nz));
  // CO
  // VERY IMPORTANT, as we need to reduce the symmetry of the derivative structure
  // getting the basis serves as a validation for the symmetry operator
  // get_full_basis = get_full_basis || derivative_structure;

  // CO START
  //  symmetry stuff
  b.ClearSymmetry(); // clear first
  // CO END

  // DX20190319 - added option to expand strictly by uniform supercell matrix - START
  // CO20190409 - note, force_supercell_matrix is PURELY for speed up purposes, the other approach should yield the same results (just slower)
  if (force_supercell_matrix && aurostd::isdiagonal(supercell)) {
    dims[1] = nx;
    dims[2] = ny;
    dims[3] = nz;
  } //[CO20190520 - fixing for dims] && !derivative_structure){dim=nx;} // if here, then nx=ny=nz
  else {
    radius = RadiusSphereLattice(b.lattice);
    dims = LatticeDimensionSphere(a.lattice, radius); //[CO20190520 - EXCESSIVE]dim=max(dims)+1;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " dims=" << dims << endl;
  }
  // DX20190319 - added option to expand strictly by uniform supercell matrix - END
  //  if(LDEBUG) cerr << "DEBUG  dims=" << dims << " " << " radius=" << radius << endl;  // DEBUG

  bool match = false;

  // CO20190409 - this pc2scMap issue is more complicated...
  // this is how we resolve: force_strict_pc2scMap
  // if force_strict_pc2scMap==false (default), then the pc2scMap returns the first of the equivalent atoms
  // this is good because of how the algorithm enumerates equivalent atoms:
  // if you know the supercell size (POCC N_HNF), then the equivalent atoms are next n_hnf atoms
  // however, for APL, we want pc2scMap to return the atom of the original primitive cell (not an equivalent atom)
  // therefore, use force_strict_pc2scMap=true
  // for force_strict_pc2scMap==true, the supercell matrix should be diagonal, otherwise you are not guaranteed to get the i==0 && j==0 && k==0 atom
  // if force_strict_pc2scMap==true and it does not find the i==0 && j==0 && k==0, it throws an error
  // CO20190114 - pc2scMap only makes sense for true supercell expansions
  // there are cases where we use GetSuperCell to convert between representations (see aflow_lattice.cpp)
  // in this case, a mapping is not really possible/useful, so clear out pc2scMap
  bool ignore_pcmap = false;
  bool pcmap = false;

  // CO20181226 - do a check if iatoms really calculated for this cell
  uint atoms_size_check = 0;
  for (size_t i = 0; i < a.iatoms.size(); i++) {
    atoms_size_check += a.iatoms[i].size();
  }
  if (atoms_size_check != a.atoms.size()) {
    a.ClearSymmetry();
  }
  // it's possible that the cell was transformed somewhere in aflow, and the symmetry was not cleared

  // CO START
  if (a.iatoms_calculated) {
    // CO START
    // create bins
    for (size_t ii = 0; ii < a.iatoms.size(); ii++) {
      b.iatoms.emplace_back(0);
    }
    for (size_t ia = 0; ia < a.iatoms.size(); ia++) {
      for (size_t iia = 0; iia < a.iatoms[ia].size(); iia++) {
        pcmap = false;
        //[CO20190520 - EXCESSIVE]for(i=-dim;i<=dim;i++) {  //[CO20200106 - close bracket for indenting]}
        //[CO20190520 - EXCESSIVE]  for(j=-dim;j<=dim;j++) {  //[CO20200106 - close bracket for indenting]}
        //[CO20190520 - EXCESSIVE]    for(k=-dim;k<=dim;k++) {  //[CO20200106 - close bracket for indenting]}
        for (i = -dims[1]; i <= dims[1]; i++) {
          for (j = -dims[2]; j <= dims[2]; j++) {
            for (k = -dims[3]; k <= dims[3]; k++) {
              atom = a.atoms[a.iatoms[ia][iia]];
              // cerr << "atom " << a.iatoms[ia][iia] << " fpos_UNrot " << atom.fpos << endl;
              cshift = ((double) i) * a.lattice(1) + ((double) j) * a.lattice(2) + ((double) k) * a.lattice(3);
              atom.cpos = atom.cpos + cshift;
              atom.fpos = b.c2f * atom.cpos; // C2F(b.lattice,atom.cpos);               // put in fractional of new basis
              //  atom.fpos=roundoff(atom.fpos);
              // cerr << "atom " << a.iatoms[ia][iia] << " fpos_rot   " << atom.fpos << endl;
              // cerr << "atom_fpos: " << atom.fpos << endl;
              // DX20180726 - if(inCell(atom.fpos)) //hard cut off
              if (uniqueAtomInCell(atom, b.atoms)) // soft cut off; then check images later //DX20180726
              { // CO20200106 - patching for auto-indenting
                // if(atom.fpos(1)>=-_eps_scell_ && atom.fpos(1)<1.0-_eps_scell_ &&
                //   atom.fpos(2)>=-_eps_scell_ && atom.fpos(2)<1.0-_eps_scell_ &&
                //   atom.fpos(3)>=-_eps_scell_ && atom.fpos(3)<1.0-_eps_scell_)      // found something inside
                //  atom=BringInCell(atom,b.lattice);
                b.num_each_type[atom.type]++; // CONVASP_MODE
                b.comp_each_type[atom.type] += atom.partial_occupation_value; // CONVASP_MODE

                // we found a new iatom
                if (b.iatoms[ia].empty()) {
                  atom.equivalent = b.atoms.size(); // reference self
                  atom.is_inequivalent = true; // iatom
                } else {
                  // eq atom
                  atom.equivalent = b.iatoms[ia][0]; // reference first iatom
                  atom.is_inequivalent = false; // equivalent atom
                }

                // ijk
                atom.ijk(1) = i;
                atom.ijk(2) = j;
                atom.ijk(3) = k;

                // DX20180726 - bring close to origin
                atom = BringCloseToOrigin(atom, b.f2c); // updates fpos/cpos/ijk

                b.atoms.push_back(atom);
                // do NOT use AddAtom(), AddAtom() rearranges per species and we need to know the mapping
                b.iatoms[ia].push_back(b.atoms.size() - 1);
                // save cshifts for fgroups later...
                match = false;
                for (size_t cf = 0; cf < cshifts.size() && !match; cf++) {
                  if (aurostd::isequal(cshift, cshifts[cf], _ZERO_TOL_)) {
                    match = true;
                  }
                }
                if (!match) {
                  cshifts.push_back(cshift);
                }

                // mapping
                sc2pcMap.push_back(a.iatoms[ia][iia]);
                // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
                if (!ignore_pcmap && !pcmap && !force_strict_pc2scMap) {
                  // if(force_strict_pc2scMap==true){  //only if i==0 && j==0 && k==0 atom
                  //   if(i==0 && j==0 && k==0){pc2scMap.push_back(b.atoms.size()-1);pcmap=true;}
                  // } else {pc2scMap.push_back(b.atoms.size()-1);pcmap=true;}
                  pc2scMap.push_back(b.atoms.size() - 1);
                  pcmap = true;
                }
                // matching cpos does not work!
                // matching by index (i,j,k) does not work because of inCell()
                // if(!i&&!j&&!k) pc2scMap.push_back(b.atoms.size()-1);
                // if(aurostd::identical(a.atoms[a.iatoms[ia][iia]].cpos,atom.cpos,1e-6)){
                //   pc2scMap[a.iatoms[ia][iia]]=b.atoms.size()-1;
                // }
              }
            }
          }
        }
        // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
        if (!ignore_pcmap && !pcmap && !force_strict_pc2scMap) {
          // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
          // if(force_strict_pc2scMap){
          //  message << "pc2scMap not found for atom[i=" << a.iatoms[ia][iia] << "]";
          //  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INDEX_MISMATCH_);
          //}
          ignore_pcmap = true;
          pc2scMap.clear();
        }
      }
    }

    // save the number of equivalents
    uint iequivalent = 0;
    for (size_t iat = 0; iat < b.atoms.size(); iat++) {
      if (b.atoms[iat].is_inequivalent) {
        b.atoms[iat].num_equivalents = b.iatoms[iequivalent].size();
        b.atoms[iat].index_iatoms = iequivalent;
        iequivalent++;
      }
    }

    b.iatoms_calculated = true;
    // CO END
  } else {
    for (size_t ia = 0; ia < a.atoms.size(); ia++) {
      pcmap = false;
      //[CO20190520 - EXCESSIVE]for(i=-dim;i<=dim;i++){ //[CO20200106 - close bracket for indenting]}
      //[CO20190520 - EXCESSIVE]  for(j=-dim;j<=dim;j++){ //[CO20200106 - close bracket for indenting]}
      //[CO20190520 - EXCESSIVE]    for(k=-dim;k<=dim;k++){ //[CO20200106 - close bracket for indenting]}
      for (i = -dims[1]; i <= dims[1]; i++) {
        for (j = -dims[2]; j <= dims[2]; j++) {
          for (k = -dims[3]; k <= dims[3]; k++) {
            atom = a.atoms[ia];
            // atom.cpos=atom.cpos+(((double)i)*a.lattice(1)+((double)j)*a.lattice(2)+((double)k)*a.lattice(3));
            cshift = ((double) i) * a.lattice(1) + ((double) j) * a.lattice(2) + ((double) k) * a.lattice(3);
            atom.cpos = atom.cpos + cshift;
            atom.fpos = b.c2f * atom.cpos; // C2F(b.lattice,atom.cpos);               // put in fractional of new basis
            //  atom.fpos=roundoff(atom.fpos);
            // DX20180726 - if(inCell(atom.fpos)) //hard cut off
            if (uniqueAtomInCell(atom, b.atoms)) // soft cut off; then check images later //DX20180726
            // if(atom.fpos(1)>=-_eps_scell_ && atom.fpos(1)<1.0-_eps_scell_ &&
            //   atom.fpos(2)>=-_eps_scell_ && atom.fpos(2)<1.0-_eps_scell_ &&
            //   atom.fpos(3)>=-_eps_scell_ && atom.fpos(3)<1.0-_eps_scell_)      // found something inside
            //  atom=BringInCell(atom,b.lattice);
            { // CO20200106 - patching for auto-indenting
              b.num_each_type[atom.type]++; // CONVASP_MODE
              b.comp_each_type[atom.type] += atom.partial_occupation_value; // CONVASP_MODE
              // ijk
              atom.ijk(1) = i;
              atom.ijk(2) = j;
              atom.ijk(3) = k;
              // DX20180726 - bring close to origin
              atom = BringCloseToOrigin(atom, b.f2c); // updates fpos/cpos/ijk
              b.atoms.push_back(atom);
              // do NOT use AddAtom(), AddAtom() rearranges per species and we need to know the mapping
              // save cshifts for fgroups later...
              match = false;
              for (size_t cf = 0; cf < cshifts.size() && !match; cf++) {
                if (aurostd::isequal(cshift, cshifts[cf], _ZERO_TOL_)) {
                  match = true;
                }
              }
              if (!match) {
                cshifts.push_back(cshift);
              }

              // mapping
              sc2pcMap.push_back(ia);
              // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
              if (!ignore_pcmap && !pcmap && !force_strict_pc2scMap) {
                // if(force_strict_pc2scMap==true){  //only if i==0 && j==0 && k==0 atom
                //   if(i==0 && j==0 && k==0){pc2scMap.push_back(b.atoms.size()-1);pcmap=true;}
                // } else {pc2scMap.push_back(b.atoms.size()-1);pcmap=true;}
                pc2scMap.push_back(b.atoms.size() - 1);
                pcmap = true;
              }
              // matching cpos does not work!
              // matching by index (i,j,k) does not work because of inCell()
              // if(!i&&!j&&!k) pc2scMap.push_back(b.atoms.size()-1);
              // if(aurostd::identical(a.atoms[ia].cpos,atom.cpos,1e-6)){
              //   pc2scMap[a.iatoms[ia][iia]]=b.atoms.size()-1;
              // }
            }
          }
        }
      }
      // ME20210506 - Strict mapping is done outside to account for non-diagonal supercells
      if (!ignore_pcmap && !pcmap && !force_strict_pc2scMap) {
        // if(force_strict_pc2scMap){
        //   message << "pc2scMap not found for atom[i=" << ia << "]";
        //   throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _INDEX_MISMATCH_);
        // }
        ignore_pcmap = true;
        pc2scMap.clear();
      }
    }
  }
  // CO END
  //  ME20210506 - The old method for force_strict_pc2scMap only works for diagonal
  //  supercells. This method is brute-force but should work for most non-diagonal
  //  cells.
  if (!ignore_pcmap && force_strict_pc2scMap) {
    pc2scMap.clear();
    const size_t pcatoms = a.atoms.size();
    pc2scMap.resize(pcatoms);
    const size_t scatoms = b.atoms.size();
    const size_t nshifts = cshifts.size();
    uint s = 0;
    uint ipc = 0;
    uint isc = 0;
    for (s = 0; s < nshifts; s++) {
      for (ipc = 0; ipc < pcatoms; ipc++) {
        for (isc = 0; isc < scatoms; isc++) {
          if (aurostd::identical(a.atoms[ipc].cpos + cshifts[s], b.atoms[isc].cpos, _FLOAT_TOL_)) {
            pc2scMap[ipc] = isc;
            break;
          }
        }
        if (isc == scatoms) {
          break;
        }
      }
      if (ipc == pcatoms) {
        break;
      }
    }
    if (s == nshifts) {
      message << "pc2scMap not found";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
  }

  b.GetStoich(); // CO20170724
  b.MakeBasis(); // need to update NUMBER and BASIS

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " sc2pcMap=" << aurostd::joinWDelimiter(sc2pcMap, " ") << endl;
    cerr << __AFLOW_FUNC__ << " pc2scMap=" << aurostd::joinWDelimiter(pc2scMap, " ") << endl;
  }

  // some check
  const double fraction = (GetVol(b.lattice) / GetVol(a.lattice));
  const double density_a = a.scale * a.scale * a.scale * std::abs(det(a.lattice)) / a.atoms.size();
  const double density_b = b.scale * b.scale * b.scale * std::abs(det(b.lattice)) / b.atoms.size();
  if (std::abs(b.atoms.size() - a.atoms.size() * fraction) > 0.1 || std::abs(density_a - density_b) > 0.001) {
    message << "Supercell has the wrong number of atoms" << endl;
    message << "b.atoms.size()     = " << b.atoms.size() << endl;
    message << "a.atoms.size()     = " << a.atoms.size() << endl;
    message << "b.scale            = " << b.scale << endl;
    message << "a.scale            = " << a.scale << endl;
    message << "fraction           = " << fraction << endl;
    message << "supercell atoms    = " << fraction * a.atoms.size() << endl;
    message << "b.density          = " << density_b << endl;
    message << "a.density          = " << density_a << endl;
    message << "b.lattice          = " << endl;
    message << b.lattice << endl;
    message << "a.lattice          = " << endl;
    message << a.lattice << endl;

    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  const bool pretend_uniform = false; // true;  //CO TEST, REMOVE ME
  get_full_basis = !pretend_uniform && get_full_basis;

  // CO START
  // for now, we focus on pgroup, fgroup, iatoms, and site symmetry (agroup), add as you need
  if (get_symmetry) {
    ofstream FileMESSAGE;
    _aflags aflags;
    _kflags kflags;
    const bool _write_ = false; // CO no verbose, annoys JJPR
    const bool osswrite = false; // CO no verbose, annoys JJPR
    bool same_pgroups = true;
    bool calculated_pgroups = false;
    const bool CALCULATE_FULL_SYMMETRY_ROBUSTLY = false; // for testing APL
    bool KRUN = true; // FORCE FULL CALC//false; //CO20181226
    // ostream& oss=cout;  //defined in macro at top of file
    if (derivative_structure) {
      KRUN = KRUN && SYM::CalculatePointGroup(FileMESSAGE, b, aflags, _write_, osswrite, cout);
      if (LDEBUG) {
        if (!KRUN) {
          cerr << "Symmetry propagation FAILED with derivative structure at point group" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at point group" << endl;
        }
      }
      same_pgroups = (KRUN && SYM::PointGroupsIdentical(a.pgroup, b.pgroup, b.sym_eps, false));
      // DX20171207 - added is_same_lattice
      calculated_pgroups = KRUN;
    }
    if ((KRUN && !pretend_uniform && !same_pgroups) || CALCULATE_FULL_SYMMETRY_ROBUSTLY) {
      if (CALCULATE_FULL_SYMMETRY_ROBUSTLY) {
        cerr << "Calculating symmetry of supercell robustly" << endl;
      }
      // we calculate earlier to see if there's a mismatch
      // SYM::CalculatePointGroup(FileMESSAGE,b,aflags,_write_,osswrite,oss);
      // do all at same sym_eps as primitive cell
      KRUN = KRUN && SYM::CalculateFactorGroup(FileMESSAGE, b, aflags, _write_, osswrite, cout);
      if (LDEBUG) {
        if (!KRUN) {
          cerr << "Symmetry propagation FAILED with derivative structure at factor group" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at factor group" << endl;
        }
      }
      KRUN = KRUN && SYM::CalculatePointGroupCrystal(FileMESSAGE, b, aflags, _write_, osswrite, cout);
      if (LDEBUG) {
        if (!KRUN) {
          cerr << "Symmetry propagation FAILED with derivative structure at point group crystal" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at point group crystal" << endl;
        }
      }
      // if(!a.iatoms_calculated){
      KRUN = KRUN && SYM::CalculateInequivalentAtoms(FileMESSAGE, b, aflags, _write_, osswrite, cout);
      // 100% necessary, new pgroups means different symmetry, different iatoms
      if (LDEBUG) {
        if (!KRUN) {
          cerr << "Symmetry propagation FAILED with derivative structure at iatoms" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at iatoms" << endl;
        }
      }
      int agroup_calculation_mode = (get_full_basis ? 0 : 1);
      if (CALCULATE_FULL_SYMMETRY_ROBUSTLY) {
        agroup_calculation_mode = 2;
      }
      //}
      // AGAIN, many fgroups, but not many pgroups for derivative structures, let's see if this is faster...
      // if(!SYM::CalculateSitePointGroup(FileMESSAGE,b,true,aflags,_write_,osswrite,oss)){  //iatoms only
      KRUN = KRUN && SYM::CalculateSitePointGroup(FileMESSAGE, b, agroup_calculation_mode, aflags, _write_, osswrite, cout);
      // don't waste time calculating basis_map for eatoms, really never use them anyway
      if (LDEBUG) {
        if (!KRUN) {
          cerr << "Symmetry propagation FAILED with derivative structure at agroup" << endl;
        } else {
          cerr << "Symmetry propagation PASSED with derivative structure at agroup" << endl;
        }
      }
      //}
      // validate that we have good symmetry here
    } else if (KRUN) {
      //////////////////////////////////////////////////////////////////////////
      // PGROUP
      if (KRUN && a.pgroup_calculated && !calculated_pgroups) {
        for (size_t i = 0; i < a.pgroup.size() && KRUN; i++) {
          pSymOp = a.pgroup[i];
          pSymOp.Uf = b.c2f * pSymOp.Uc * b.f2c;
          pSymOp.basis_atoms_map.clear();
          pSymOp.basis_types_map.clear();
          pSymOp.basis_map_calculated = false;
          // CO+DX, getFullSymBasis does not make sense for just rotations
          // if(get_full_basis){
          //   KRUN = KRUN && SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,pSymOp,false,skew,b.sym_eps,pSymOp.basis_atoms_map,pSymOp.basis_types_map);
          //   if(LDEBUG) {
          //     if(!KRUN){
          //       cerr << "Symmetry propagation FAILED with uniform supercell structure at point group" << endl;
          //     } else {
          //       cerr << "Symmetry propagation PASSED with uniform supercell structure at point group" << endl;
          // }
          //   }
          //   pSymOp.basis_map_calculated=KRUN;
          // }
          // if(KRUN){b.pgroup.push_back(pSymOp);}
          if (KRUN) {
            SYM::AddSymmetryToStructure(b, pSymOp.Uc, pSymOp.Uf, pSymOp.ctau, pSymOp.ftau, pSymOp.ctrasl, pSymOp.ftrasl, pSymOp.basis_atoms_map, pSymOp.basis_types_map, pSymOp.basis_map_calculated, _PGROUP_, false);
          } // CO20170706 - make sure quaternion is updated
        }
        // if(KRUN && b.pgroup.size()){
        // b.pgroup_calculated=true;
        b.pgroup_calculated = (KRUN && !b.pgroup.empty());
        //}
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // FGROUP
      if (KRUN && a.fgroup_calculated) {
        for (size_t fg = 0; fg < a.fgroup.size() && KRUN; fg++) {
          for (size_t cs = 0; cs < cshifts.size() && KRUN; cs++) {
            fSymOp = a.fgroup[fg];
            fSymOp.ctau = fSymOp.ctau + cshifts[cs];
            fSymOp.ftau = b.c2f * fSymOp.ctau;
            if (inCell(fSymOp.ftau)) { // DX CHANGE HERE; NO MORE TOL_ABC_RES
              // We have to correct the Uf for each symop since we have changed the lattice...
              fSymOp.Uf = b.c2f * fSymOp.Uc * b.f2c;
              fSymOp.basis_atoms_map.clear();
              fSymOp.basis_types_map.clear();
              fSymOp.basis_map_calculated = false;
              for (size_t iii = 0; iii < b.atoms.size(); iii++) {
                fSymOp.basis_atoms_map.push_back(0);
                fSymOp.basis_types_map.push_back(0);
              }
              fSymOp.basis_map_calculated = false;
              // calculate basis_atoms_map and basis_types_map
              if (get_full_basis) {
                // NOPE, we will calculate
                KRUN = KRUN && SYM::getFullSymBasis(b.atoms, b.lattice, b.c2f, b.f2c, fSymOp, true, skew, b.sym_eps, fSymOp.basis_atoms_map, fSymOp.basis_types_map);
                if (LDEBUG) {
                  if (!KRUN) {
                    cerr << "Symmetry propagation FAILED with uniform supercell structure at factor group" << endl;
                  } else {
                    cerr << "Symmetry propagation PASSED with uniform supercell structure at factor group" << endl;
                  }
                }
                // if(!SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,fSymOp,true,skew,b.sym_eps,fSymOp.basis_atoms_map,fSymOp.basis_types_map)) {
                // cerr << "Unable to find atom/types basis for fgroup" << endl;
                // cerr << fSymOp << endl;
                // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
                // KRUN = false;
                // }
                fSymOp.basis_map_calculated = KRUN;
              }
              // if(KRUN){b.fgroup.push_back(fSymOp);}
              if (KRUN) {
                SYM::AddSymmetryToStructure(b, fSymOp.Uc, fSymOp.Uf, fSymOp.ctau, fSymOp.ftau, fSymOp.ctrasl, fSymOp.ftrasl, fSymOp.basis_atoms_map, fSymOp.basis_types_map, fSymOp.basis_map_calculated, _FGROUP_, false);
              } // CO20170706 - make sure quaternion is updated
            }
          }
        }
        // if(b.fgroup.size()){
        //   b.fgroup_calculated=true;   //there are more, but these are the important ones
        b.fgroup_calculated = (KRUN && !b.fgroup.empty()); // there are more, but these are the important ones
        //}
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // PGROUP_XTAL
      if (KRUN && a.pgroup_xtal_calculated) {
        for (size_t i = 0; i < a.pgroup_xtal.size() && KRUN; i++) {
          pSymOp = a.pgroup_xtal[i];
          pSymOp.Uf = b.c2f * pSymOp.Uc * b.f2c;
          pSymOp.basis_atoms_map.clear();
          pSymOp.basis_types_map.clear();
          pSymOp.basis_map_calculated = false;
          // CO+DX, getFullSymBasis does not make sense for just rotations
          // if(get_full_basis){
          //   KRUN = KRUN && SYM::getFullSymBasis(b.atoms,b.lattice,b.c2f,b.f2c,pSymOp,false,skew,b.sym_eps,pSymOp.basis_atoms_map,pSymOp.basis_types_map);
          //   if(LDEBUG) {
          //     if(!KRUN){
          //       cerr << "Symmetry propagation FAILED with uniform supercell structure at point group crystal" << endl;
          //     } else {
          //       cerr << "Symmetry propagation PASSED with uniform supercell structure at point group crystal" << endl;
          //     }
          //   }
          //   pSymOp.basis_map_calculated=KRUN;
          // }
          // if(KRUN){b.pgroup_xtal.push_back(pSymOp);}
          if (KRUN) {
            SYM::AddSymmetryToStructure(b, pSymOp.Uc, pSymOp.Uf, pSymOp.ctau, pSymOp.ftau, pSymOp.ctrasl, pSymOp.ftrasl, pSymOp.basis_atoms_map, pSymOp.basis_types_map, pSymOp.basis_map_calculated, _PGROUP_XTAL_, false);
          } // CO20170706 - make sure quaternion is updated
        }
        // if(KRUN && b.pgroup_xtal.size()){
        // b.pgroup_xtal_calculated=true;
        b.pgroup_xtal_calculated = (KRUN && !b.pgroup_xtal.empty());
        //}
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // IF A.IATOMS, then we already propagated, otherwise, just run full routine for iatoms
      if (KRUN && !a.iatoms_calculated) {
        KRUN = KRUN && SYM::CalculateInequivalentAtoms(FileMESSAGE, b, aflags, _write_, osswrite, cout); // do all atoms
        if (LDEBUG) {
          if (!KRUN) {
            cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
          } else {
            cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
          }
        }
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // AGROUP
      if (KRUN && a.agroup_calculated) {
        if (b.iatoms_calculated && b.fgroup_calculated && get_full_basis) {
          // the only way we get a speed up is if we can use fgroup.basis_atoms_map
          // xstructure bb=b;
          deque<_atom> b_atoms = b.atoms;
          xvector<double> origin(3);
          xvector<double> frigin(3);
          // create space for agroups
          for (size_t iia = 0; iia < b.atoms.size(); iia++) {
            b.agroup.emplace_back(0);
          }
          for (size_t ia = 0; ia < b.iatoms.size() && KRUN; ia++) {
            iat = b.iatoms[ia][0];
            // let's recycle what we have to improve speed
            origin = b.atoms[iat].cpos;
            frigin = b.atoms[iat].fpos;
            for (size_t ii = 0; ii < b.atoms.size(); ii++) {
              // go back to original first, then subtract new origin
              b_atoms[ii].cpos = b.atoms[ii].cpos - origin;
              b_atoms[ii].fpos = b.atoms[ii].fpos - frigin;
              // now bring in cell
              b_atoms[ii] = BringInCell(b_atoms[ii], b.lattice);
            }
            // bb.ShiftOriginToAtom(iat);
            // bb.BringInCell();

            // IATOMS ONLY
            for (size_t iia = 0; iia < a.agroup[sc2pcMap[iat]].size() && KRUN; iia++) {
              aSymOp = a.agroup[sc2pcMap[iat]][iia];
              // We have to correct the Uf for each symop since we have changed the lattice...
              aSymOp.Uf = b.c2f * aSymOp.Uc * b.f2c;
              // no longer necessary since we force a basis map calculation
              // aSymOp.basis_atoms_map.clear();
              // aSymOp.basis_types_map.clear();
              // aSymOp.basis_map_calculated=false;
              // for(size_t iii=0;iii<b.atoms.size();iii++){
              // aSymOp.basis_atoms_map.push_back(0);
              // aSymOp.basis_types_map.push_back(0);
              // }
              // fSymOp.basis_map_calculated=false;
              KRUN = KRUN && SYM::getFullSymBasis(b_atoms, b.lattice, b.c2f, b.f2c, aSymOp, true, skew, b.sym_eps, aSymOp.basis_atoms_map, aSymOp.basis_types_map);
              if (LDEBUG) {
                if (!KRUN) {
                  cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
                } else {
                  cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
                }
              }
              // if(!SYM::getFullSymBasis(b_atoms,b.lattice,b.c2f,b.f2c,aSymOp,true,skew,b.sym_eps,aSymOp.basis_atoms_map,aSymOp.basis_types_map)) {
              // cerr << "Unable to find atom/types basis for agroup" << endl;
              // cerr << aSymOp << endl;
              // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
              // KRUN = false;
              // }
              aSymOp.basis_map_calculated = KRUN;
              // if(KRUN){b.agroup[iat].push_back(aSymOp);}
              if (KRUN) {
                SYM::AddSymmetryToStructure(b, iat, aSymOp.Uc, aSymOp.Uf, aSymOp.ctau, aSymOp.ftau, aSymOp.ctrasl, aSymOp.ftrasl, aSymOp.basis_atoms_map, aSymOp.basis_types_map, aSymOp.basis_map_calculated, _AGROUP_, false);
              } // CO20170706 - make sure quaternion is updated
            }
          }

          // EATOMS FOLLOW
          KRUN = KRUN && SYM::CalculateSitePointGroup_EquivalentSites(b, get_full_basis, b.sym_eps);
          if (LDEBUG) {
            if (!KRUN) {
              cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup equivalent" << endl;
            } else {
              cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup equivalent" << endl;
            }
          }
          // KRUN = false;
          // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unable to propagate site symmetry to equivalent atoms.",_RUNTIME_ERROR_);
        } else {
          // can be faster than procedure above because there are MANY fgroups
          KRUN = KRUN && SYM::CalculateSitePointGroup(FileMESSAGE, b, 1, aflags, _write_, osswrite, cout);
          // we already know get_full_basis==false, so don't waste time calculating for eatoms
          if (LDEBUG) {
            if (!KRUN) {
              cerr << "Symmetry propagation FAILED with uniform supercell structure at agroup" << endl;
            } else {
              cerr << "Symmetry propagation PASSED with uniform supercell structure at agroup" << endl;
            }
          }
        }
        // if(b.agroup.size() && b.agroup[0].size()){      //just a fast check to see we have agroups somewhere (we should always get identity)
        b.agroup_calculated = (KRUN && !b.agroup[0].empty());
        // just a fast check to see we have agroups somewhere (we should always get identity)
        //   b.agroup_calculated=true;
        // }
      }
      //////////////////////////////////////////////////////////////////////////
    }
    if (!KRUN) {
      cout << (aflags.QUIET ? "" : "00000  MESSAGE ") << "SUPERCELL Symmetry propagation FAILED" << Message(__AFLOW_FILE__, aflags) << endl;
      cout << (aflags.QUIET ? "" : "00000  MESSAGE ") << "SUPERCELL Symmetry retrying with symmetry scan" << Message(__AFLOW_FILE__, aflags) << endl;
      b.ClearSymmetry(); // CO20181226
      pflow::PerformFullSymmetry(b, FileMESSAGE, aflags, kflags, osswrite, cout);
      // FOOLPROOF!!!!!!!!!
      // no need for krun here with force_perform: if it fails with full scan, it will calculate at default tolerance and keep going!
    }
  }

  // cerr << "prim: " << a.pgroup.size() << ", sup: " << b.pgroup.size() << endl;    //CO REMOVE

  // CO END

  b.ReScale(aa.scale); // the nuclear option, the only way not to mess around with scale EVERYWHERE
  return b;
}

// CO END

xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell, vector<int>& sc2pcMap, vector<int>& pc2scMap, bool get_symmetry, bool get_full_basis, bool force_supercell_matrix, bool force_strict_pc2scMap)
// DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
// xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell)
{ // CO20200106 - patching for auto-indenting
  const xmatrix<double> _supercell(3, 3);
  if (supercell.rows == 9) {
    _supercell(1, 1) = supercell(1);
    _supercell(1, 2) = supercell(2);
    _supercell(1, 3) = supercell(3);
    _supercell(2, 1) = supercell(4);
    _supercell(2, 2) = supercell(5);
    _supercell(2, 3) = supercell(6);
    _supercell(3, 1) = supercell(7);
    _supercell(3, 2) = supercell(8);
    _supercell(3, 3) = supercell(9);
    return GetSuperCell(a, _supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
    // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  }
  if (supercell.rows == 3) {
    _supercell(1, 1) = supercell(1);
    _supercell(2, 2) = supercell(2);
    _supercell(3, 3) = supercell(3);
    return GetSuperCell(a, _supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
    // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  }
  stringstream message;
  message << "Matrix must have 9 or 3 elements";
  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
}

xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell, vector<int>& sc2pcMap, vector<int>& pc2scMap, bool get_symmetry, bool get_full_basis, bool force_supercell_matrix, bool force_strict_pc2scMap)
// DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
// xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell)
{ // CO20200106 - patching for auto-indenting
  const xmatrix<double> _supercell(3, 3);
  if (supercell.rows == 9) {
    _supercell(1, 1) = supercell(1);
    _supercell(1, 2) = supercell(2);
    _supercell(1, 3) = supercell(3);
    _supercell(2, 1) = supercell(4);
    _supercell(2, 2) = supercell(5);
    _supercell(2, 3) = supercell(6);
    _supercell(3, 1) = supercell(7);
    _supercell(3, 2) = supercell(8);
    _supercell(3, 3) = supercell(9);
    return GetSuperCell(a, _supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
    // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  }
  if (supercell.rows == 3) {
    _supercell(1, 1) = supercell(1);
    _supercell(2, 2) = supercell(2);
    _supercell(3, 3) = supercell(3);
    return GetSuperCell(a, _supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
    // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
  }
  stringstream message;
  message << "Matrix must have 9 or 3 elements";
  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
}

xstructure GetSuperCell(const xstructure& a,
                        const int& sc11,
                        const int& sc12,
                        const int& sc13,
                        const int& sc21,
                        const int& sc22,
                        const int& sc23,
                        const int& sc31,
                        const int& sc32,
                        const int& sc33,
                        vector<int>& sc2pcMap,
                        vector<int>& pc2scMap,
                        bool get_symmetry,
                        bool get_full_basis,
                        bool force_supercell_matrix,
                        bool force_strict_pc2scMap)
// DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
// xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33)
{ // CO20200106 - patching for auto-indenting
  xmatrix<double> _supercell(3, 3);
  _supercell.clear();
  _supercell(1, 1) = (double) sc11;
  _supercell(1, 2) = (double) sc12;
  _supercell(1, 3) = (double) sc13;
  _supercell(2, 1) = (double) sc21;
  _supercell(2, 2) = (double) sc22;
  _supercell(2, 3) = (double) sc23;
  _supercell(3, 1) = (double) sc31;
  _supercell(3, 2) = (double) sc32;
  _supercell(3, 3) = (double) sc33;
  return GetSuperCell(a, _supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
  // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const int& sc1, const int& sc2, const int& sc3, vector<int>& sc2pcMap, vector<int>& pc2scMap, bool get_symmetry, bool get_full_basis, bool force_supercell_matrix, bool force_strict_pc2scMap)
// DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
// xstructure GetSuperCell(const xstructure& a, const int& sc1,const int& sc2,const int& sc3)
{ // CO20200106 - patching for auto-indenting
  xmatrix<double> _supercell(3, 3);
  _supercell.clear();
  _supercell(1, 1) = (double) sc1;
  _supercell(2, 2) = (double) sc2;
  _supercell(3, 3) = (double) sc3;
  return GetSuperCell(a, _supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
  // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

// CO START
xstructure GetSuperCell(const xstructure& aa, const xmatrix<double>& supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  const bool get_symmetry = false;
  const bool get_full_basis = false;
  const bool force_supercell_matrix = false; // DX20190319 - added force_supercell_matrix
  const bool force_strict_pc2scMap = false; // CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(aa, supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
  // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const xvector<double>& supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  const bool get_symmetry = false;
  const bool get_full_basis = false;
  const bool force_supercell_matrix = false; // DX20190319 - added force_supercell_matrix
  const bool force_strict_pc2scMap = false; // CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(a, supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
  // DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const xvector<int>& supercell) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  const bool get_symmetry = false;
  const bool get_full_basis = false;
  const bool force_supercell_matrix = false; // DX20190319 - added force_supercell_matrix
  const bool force_strict_pc2scMap = false; // CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(a, supercell, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
  // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const int& sc11, const int& sc12, const int& sc13, const int& sc21, const int& sc22, const int& sc23, const int& sc31, const int& sc32, const int& sc33) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  const bool get_symmetry = false;
  const bool get_full_basis = false;
  const bool force_supercell_matrix = false; // DX20190319 - added force_supercell_matrix
  const bool force_strict_pc2scMap = false; // CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(a, sc11, sc12, sc13, sc21, sc22, sc23, sc31, sc32, sc33, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
  // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
}

xstructure GetSuperCell(const xstructure& a, const int& sc1, const int& sc2, const int& sc3) {
  vector<int> sc2pcMap;
  vector<int> pc2scMap;
  const bool get_symmetry = false;
  const bool get_full_basis = false;
  const bool force_supercell_matrix = false; // DX20190319 - added force_supercell_matrix
  const bool force_strict_pc2scMap = false; // CO20190409 - added force_strict_pc2scMap
  return GetSuperCell(a, sc1, sc2, sc3, sc2pcMap, pc2scMap, get_symmetry, get_full_basis, force_supercell_matrix, force_strict_pc2scMap);
  // DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
}

// CO END

// ***************************************************************************
// Function ClearSymmetry
// ***************************************************************************
void xstructure::ClearSymmetry() {
  const bool LDEBUG = (false || XHOST.DEBUG); // DX20210406

  // PGROUP ----------------------------
  pgroup.clear(); // just initialize
  pgroup_calculated = false;
  bravais_lattice_lattice_type = "";
  bravais_lattice_lattice_variation_type = "";
  bravais_lattice_lattice_system = ""; // DX20210430 - missing
  // PGROUP_XTAL ----------------------------
  pgroup_xtal.clear(); // just initialize
  pgroup_xtal_calculated = false;
  crystal_family = "";
  crystal_system = "";
  point_group_crystal_class = "";
  point_group_Shoenflies = "";
  point_group_Hermann_Mauguin = "";
  point_group_orbifold = "";
  point_group_type = "";
  point_group_order = "";
  point_group_structure = "";
  pearson_symbol = ""; // DX20210430 - missing
  bravais_lattice_type = "";
  bravais_lattice_variation_type = "";
  bravais_lattice_system = ""; // DX20210430 - missing
  // PGROUPK_PATTERSON ---------------------------- //DX20200129
  pgroupk_Patterson.clear(); // just initialize
  pgroupk_Patterson_calculated = false;
  // PGROUPK ----------------------------
  pgroupk.clear(); // just initialize
  pgroupk_calculated = false;
  reciprocal_lattice_type = "";
  reciprocal_lattice_variation_type = ""; // DX20210430 - missing
  // PGROUPK_XTAL ----------------------------
  pgroupk_xtal.clear(); // just initialize //DX20171205 - Added pgroupk_xtal
  pgroupk_xtal_calculated = false; // DX20171205 - Added pgroupk_xtal
  // FGROUP ----------------------------
  fgroup.clear(); // just initialize
  fgroup_calculated = false;
  // SGROUP ----------------------------
  sgroup_radius = -_calculate_symmetry_default_sgroup_radius_; // symmetry not calculated
  sgroup_radius_dims.clear();
  sgroup.clear(); // just initialize
  sgroup_calculated = false;
  // SITE POINT GROUP ------------------
  agroup_calculated = false;
  for (size_t i = 0; i < agroup.size(); i++) {
    agroup[i].clear();
  }
  agroup.clear();
  // INEQUIVALENT ATOMS ----------------
  iatoms_calculated = false;
  for (size_t i = 0; i < iatoms.size(); i++) {
    iatoms[i].clear();
  }
  iatoms.clear();
  for (size_t i = 0; i < atoms.size(); i++) {
    atoms[i].ClearSymmetry();
  } // CO20190219
  // SUPERLATTICE //DX20210430 - missing
  bravais_superlattice_lattice.clear();
  bravais_superlattice_type = "";
  bravais_superlattice_variation_type = "";
  bravais_superlattice_system = "";
  pearson_symbol_superlattice = "";

  if (LDEBUG) {
    cerr << XPID << "xstructure::ClearSymmetry(): All symmetry attributes have been cleared." << endl;
  } // DX20210406
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//       It contains consistency checks for the symmetry analysis.
//  ***************************************************************************
//  Function CalculateSymmetry
//  ***************************************************************************
bool xstructure::CalculateSymmetry(bool ossverbose, double radius) {
  ofstream FileDevNull("/dev/null");

  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  _kflags kflags;
  pflow::defaultKFlags4SymWrite(kflags, false);
  pflow::defaultKFlags4SymCalc(kflags, true);

  (*this).LatticeReduction_avoid = true;
  (*this).sgroup_radius = radius;

  return pflow::PerformFullSymmetry(*this, FileDevNull, aflags, kflags, ossverbose, cout);

  // SYM::CalculatePointGroup(FileDevNull,*this,aflags,false,ossverbose,cout);
  // SYM::CalculatePointGroupKLattice(FileDevNull,*this,aflags,false,ossverbose,cout);
  // SYM::CalculateSitePointGroup(FileDevNull,*this,aflags,false,ossverbose,cout);
  // SYM::CalculateFactorGroup(FileDevNull,*this,aflags,false,ossverbose,cout);
  // SYM::CalculateInequivalentAtoms(FileDevNull,*this,aflags,false,ossverbose,cout);
  // SYM::CalculateSpaceGroup(FileDevNull,*this,aflags,false,ossverbose,cout);
}

bool xstructure::CalculateSymmetry() {
  LatticeReduction_avoid = true;
  return CalculateSymmetry(false, _calculate_symmetry_default_sgroup_radius_ * MaxStructureLattice(*this));
}

bool CalculateSymmetry(xstructure& str, bool ossverbose, ostream& cout, bool fffverbose, double radius) {
  ofstream FileDevNull("/dev/null");

  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  _kflags kflags;
  pflow::defaultKFlags4SymWrite(kflags, fffverbose);
  pflow::defaultKFlags4SymCalc(kflags, true);

  str.LatticeReduction_avoid = true;
  str.sgroup_radius = radius;

  return pflow::PerformFullSymmetry(str, FileDevNull, aflags, kflags, ossverbose, cout);
  // SYM::CalculatePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  // SYM::CalculatePointGroupKLattice(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  // SYM::CalculateSitePointGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  // SYM::CalculateFactorGroup(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  // SYM::CalculateInequivalentAtoms(FileDevNull,str,aflags,fffverbose,ossverbose,oss);
  // SYM::CalculateSpaceGroup(FileDevNull,str,aflags,true,ossverbose,oss);
}

bool CalculateSymmetry(xstructure& str, bool ossverbose, ostream& cout, bool fffverbose) {
  str.LatticeReduction_avoid = true;
  return CalculateSymmetry(str, ossverbose, cout, fffverbose, _calculate_symmetry_default_sgroup_radius_ * MaxStructureLattice(str));
}

bool CalculateSymmetry(xstructure& str, bool ossverbose, ostream& cout, double radius) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  return CalculateSymmetry(str, ossverbose, cout, false, radius);
}

bool CalculateSymmetry(xstructure& str, bool ossverbose, double radius) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  return CalculateSymmetry(str, ossverbose, cout, false, radius);
}

bool CalculateSymmetry(xstructure& str, double radius) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  return CalculateSymmetry(str, false, cout, false, radius);
}

bool CalculateSymmetry(xstructure& str, bool ossverbose) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  return CalculateSymmetry(str, ossverbose, cout, false, _calculate_symmetry_default_sgroup_radius_ * MaxStructureLattice(str));
}

bool CalculateSymmetry(xstructure& str) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  return CalculateSymmetry(str, false, cout, false, _calculate_symmetry_default_sgroup_radius_ * MaxStructureLattice(str));
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//       It contains consistency checks for the symmetry analysis.
//  ***************************************************************************
//  Function CalculateSymmetryPointGroup
//  ***************************************************************************
void xstructure::CalculateSymmetryPointGroup(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  LatticeReduction_avoid = true;
  SYM::CalculatePointGroup(FileDevNull, *this, aflags, false, ossverbose, cout);
}

void xstructure::CalculateSymmetryPointGroup() {
  LatticeReduction_avoid = true;
  CalculateSymmetryPointGroup(false);
}

void CalculateSymmetryPointGroup(xstructure& str, bool ossverbose, ostream& cout, bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  str.LatticeReduction_avoid = true;
  SYM::CalculatePointGroup(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
}

void CalculateSymmetryPointGroup(xstructure& str, bool ossverbose, ostream& cout) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroup(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroup(xstructure& str, bool ossverbose) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroup(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroup(xstructure& str) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroup(str, false, cout, false);
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//       It contains consistency checks for the symmetry analysis.
//  ***************************************************************************
//  Function CalculateSymmetryPointGroupCrystal
//  ***************************************************************************
void xstructure::CalculateSymmetryPointGroupCrystal(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  LatticeReduction_avoid = true;
  if (pgroup_calculated == false) {
    SYM::CalculatePointGroup(FileDevNull, *this, aflags, false, ossverbose, cout);
  }
  if (fgroup_calculated == false) {
    SYM::CalculateFactorGroup(FileDevNull, *this, aflags, false, ossverbose, cout);
  }
  SYM::CalculatePointGroupCrystal(FileDevNull, *this, aflags, false, ossverbose, cout);
}

void xstructure::CalculateSymmetryPointGroupCrystal() {
  LatticeReduction_avoid = true;
  if (pgroup_calculated == false) {
    CalculateSymmetryPointGroup(false);
  }
  if (fgroup_calculated == false) {
    CalculateSymmetryFactorGroup(false);
  }
  CalculateSymmetryPointGroupCrystal(false);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str, bool ossverbose, ostream& cout, bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  str.LatticeReduction_avoid = true;
  if (str.pgroup_calculated == false) {
    SYM::CalculatePointGroup(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
  }
  if (str.fgroup_calculated == false) {
    SYM::CalculateFactorGroup(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
  }
  SYM::CalculatePointGroupCrystal(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str, bool ossverbose, ostream& cout) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupCrystal(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str, bool ossverbose) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupCrystal(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroupCrystal(xstructure& str) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  if (str.pgroup_calculated == false) {
    CalculateSymmetryPointGroup(str, false, cout, false);
  }
  if (str.fgroup_calculated == false) {
    CalculateSymmetryFactorGroup(str, false, cout, false);
  }
  CalculateSymmetryPointGroupCrystal(str, false, cout, false);
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//       It contains consistency checks for the symmetry analysis.
//  ***************************************************************************
//  Function CalculateSymmetryFactorGroup
//  ***************************************************************************
void xstructure::CalculateSymmetryFactorGroup(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  LatticeReduction_avoid = true;
  SYM::CalculatePointGroup(FileDevNull, *this, aflags, false, ossverbose, cout);
  SYM::CalculateFactorGroup(FileDevNull, *this, aflags, false, ossverbose, cout);
}

void xstructure::CalculateSymmetryFactorGroup() {
  LatticeReduction_avoid = true;
  CalculateSymmetryFactorGroup(false);
}

void CalculateSymmetryFactorGroup(xstructure& str, bool ossverbose, ostream& cout, bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  str.LatticeReduction_avoid = true;
  SYM::CalculatePointGroup(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
  SYM::CalculateFactorGroup(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
}

void CalculateSymmetryFactorGroup(xstructure& str, bool ossverbose, ostream& cout) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryFactorGroup(str, ossverbose, cout, false);
}

void CalculateSymmetryFactorGroup(xstructure& str, bool ossverbose) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryFactorGroup(str, ossverbose, cout, false);
}

void CalculateSymmetryFactorGroup(xstructure& str) {
  str.LatticeReduction_avoid = true;
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryFactorGroup(str, false, cout, false);
}

// DX - Consider using pflow::CalculateFullSymmetry in aflow_aconvasp_main.cpp.
//       It contains consistency checks for the symmetry analysis.
// ME20200114 - made capitalization more consistent with other functions
//  ***************************************************************************
//  Function CalculateSymmetryPointGroupKLattice
//  ***************************************************************************
void xstructure::CalculateSymmetryPointGroupKLattice(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  // LatticeReduction_avoid=true;  // not necssary for klattice
  SYM::CalculatePointGroupKLattice(FileDevNull, *this, aflags, false, ossverbose, cout);
}

void xstructure::CalculateSymmetryPointGroupKLattice() {
  // LatticeReduction_avoid=true;  // not necssary for klattice
  CalculateSymmetryPointGroupKLattice(false);
}

void CalculateSymmetryPointGroupKLattice(xstructure& str, bool ossverbose, ostream& cout, bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  // str.LatticeReduction_avoid=true;  // not necssary for klattice
  SYM::CalculatePointGroupKLattice(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
}

void CalculateSymmetryPointGroupKLattice(xstructure& str, bool ossverbose, ostream& cout) {
  // str.LatticeReduction_avoid=true;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKLattice(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroupKLattice(xstructure& str, bool ossverbose) {
  //  str.LatticeReduction_avoid=true;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKLattice(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroupKLattice(xstructure& str) {
  //  str.LatticeReduction_avoid=true;  // not necssary for klattice
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKLattice(str, false, cout, false);
}

// ME20200114 - added missing function
//  ***************************************************************************
//  Function CalculateSymmetryPointGroupKCrystal
//  ***************************************************************************
void xstructure::CalculateSymmetryPointGroupKCrystal(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  SYM::CalculatePointGroupKCrystal(FileDevNull, *this, aflags, false, ossverbose, cout);
}

void xstructure::CalculateSymmetryPointGroupKCrystal() {
  CalculateSymmetryPointGroupKCrystal(false);
}

void CalculateSymmetryPointGroupKCrystal(xstructure& str, bool ossverbose, ostream& cout, bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  SYM::CalculatePointGroupKCrystal(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
}

void CalculateSymmetryPointGroupKCrystal(xstructure& str, bool ossverbose, ostream& cout) {
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKCrystal(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroupKCrystal(xstructure& str, bool ossverbose) {
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKCrystal(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroupKCrystal(xstructure& str) {
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKCrystal(str, false, cout, false);
}

// ME20200129
//  ***************************************************************************
//  Function CalculateSymmetryPointGroupKPatterson
//  ***************************************************************************

void xstructure::CalculateSymmetryPointGroupKPatterson(bool ossverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  SYM::CalculatePointGroupKPatterson(FileDevNull, *this, aflags, false, ossverbose, cout);
}

void xstructure::CalculateSymmetryPointGroupKPatterson() {
  CalculateSymmetryPointGroupKPatterson(false);
}

void CalculateSymmetryPointGroupKPatterson(xstructure& str, bool ossverbose, ostream& cout, bool fffverbose) {
  ofstream FileDevNull("/dev/null");
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  SYM::CalculatePointGroupKPatterson(FileDevNull, str, aflags, fffverbose, ossverbose, cout);
}

void CalculateSymmetryPointGroupKPatterson(xstructure& str, bool ossverbose, ostream& cout) {
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKPatterson(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroupKPatterson(xstructure& str, bool ossverbose) {
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKPatterson(str, ossverbose, cout, false);
}

void CalculateSymmetryPointGroupKPatterson(xstructure& str) {
  _aflags aflags;
  aflags.Directory = "./";
  aflags.QUIET = true;
  CalculateSymmetryPointGroupKPatterson(str, false, cout, false);
}

// ***************************************************************************
// Function fixEmptyAtomNames()
// ***************************************************************************
void xstructure::fixEmptyAtomNames(bool force_fix) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  if (species.size() == species_pp.size()) { // CO20190218
    for (size_t itype = 0; itype < species.size(); itype++) {
      if ((force_fix || species[itype].empty()) && !species_pp.at(itype).empty()) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " species_pp.at(" << itype << ")=" << species_pp.at(itype) << endl;
        }
        species[itype] = species_pp.at(itype);
        // KBIN::VASP_PseudoPotential_CleanName(species_pp.at(itype));  //CO20181226 KEEP PP INFO if available (auto aflow.in)
      }
    }
  } // cormac I`ll write a short pflow for this stuff
  if (species.size() == num_each_type.size()) {
    int iatom = 0;
    for (size_t itype = 0; itype < num_each_type.size(); itype++) {
      const string s = string(species.at(itype));
      species.at(itype) = s;
      for (int j = 0; j < num_each_type[itype]; j++) {
        atoms.at(iatom).name = s; // CONVASP_MODE
        atoms.at(iatom).CleanName();
        atoms.at(iatom).CleanSpin();
        atoms.at(iatom).name_is_given = true;
        iatom++;
      }
    }
  }
}

// ***************************************************************************
// Function buildGenericTitle()
// ***************************************************************************
void xstructure::buildGenericTitle(bool vasp_input, bool force_fix) {
  if (vasp_input) {
    pflow::fixEmptyAtomNames(*this, force_fix);
  }
  title.clear();
  title = getGenericTitleXStructure(*this, false); // no latex
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]uint iat=0;
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]//if any names missing from atoms, lets use generic names
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]bool atom_names=true;
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]for(size_t i=0;i<atoms.size()&&atom_names;i++){if(atoms[i].name.empty()){atom_names=false;}} //CO20180316 - use pp names
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]for(size_t itype=0;itype<num_each_type.size();itype++){
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]  for(uint j=0;j<(uint) num_each_type.at(itype);j++) {
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]    if(j==0){
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]      if(atom_names){title+=atoms.at(iat).name;} //CO20180316 - use pp names
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]      else {title+=char('A'+itype);}
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]      title+=aurostd::utype2string(num_each_type.at(itype));
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]    }
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]    iat++;
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]  }
  //[CO20190418 - MOVED TO aflow_pflow_functions.cpp]}
}

// ***************************************************************************
// Function platon2print
// ***************************************************************************
string xstructure::platon2print(bool P_EQUAL, bool P_EXACT, double P_ang, double P_d1, double P_d2, double P_d3) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream cout;
  cout.setf(std::ios::fixed, std::ios::floatfield);
  xstructure str = *this;

  // Deal with too small volume problem
  if (str.GetVolume() < (double) PLATON_MIN_VOLUME_PER_ATOM * str.atoms.size()) {
    str.SetVolume(PLATON_MIN_VOLUME_PER_ATOM * str.atoms.size());
    if (LDEBUG) {
      cerr << "platon2print:  PLATON FIXED VOLUME PER ATOM = " << PLATON_MIN_VOLUME_PER_ATOM << endl;
    }
  }
  if (LDEBUG) {
    cerr << "platon2print: volume=" << str.GetVolume() << endl;
  }
  // Set scale to 1 so you don't need to rescale coordinates.
  str.ReScale(1.0);
  std::vector<string> el_names(7);
  int k;
  el_names[0] = "Ag";
  el_names[1] = "Zr";
  el_names[2] = "Cd";
  el_names[3] = "Mo";
  el_names[4] = "Fe";
  el_names[5] = "W";
  el_names[6] = "O";
  // Print out data in Platon format

  cout << "TITL " << str.title << endl;
  cout.precision(8);
  cout << "CELL " << str.scale * aurostd::modulus(str.lattice(1, 1), str.lattice(1, 2), str.lattice(1, 3)) << " " << str.scale * aurostd::modulus(str.lattice(2, 1), str.lattice(2, 2), str.lattice(2, 3)) << " "
       << str.scale * aurostd::modulus(str.lattice(3, 1), str.lattice(3, 2), str.lattice(3, 3)) << " "
       << aurostd::angle(str.lattice(2, 1), str.lattice(2, 2), str.lattice(2, 3), str.lattice(3, 1), str.lattice(3, 2), str.lattice(3, 3)) << " "
       << aurostd::angle(str.lattice(1, 1), str.lattice(1, 2), str.lattice(1, 3), str.lattice(3, 1), str.lattice(3, 2), str.lattice(3, 3)) << " "
       << aurostd::angle(str.lattice(1, 1), str.lattice(1, 2), str.lattice(1, 3), str.lattice(2, 1), str.lattice(2, 2), str.lattice(2, 3)) << endl;
  // oss << str.atoms.size() << endl;
  cout.precision(10);
  k = 0;
  if (str.num_each_type.size() != 2) {
    for (size_t i = 0; i < str.num_each_type.size(); i++) {
      for (int j = 1; j <= str.num_each_type[i]; j++) {
        //      printf("%c%i ",i+66,j);
        if (str.atoms.at(k).name_is_given) {
          cout << str.atoms.at(k).cleanname << j << " ";
        } else { // Must make up a name
          if (i > (el_names.size() - 1)) { // No more default names so make everything W
            cout << "W" << j << " ";
          } else { // Use default names
            cout << el_names[i] << j << " ";
          }
        }
        cout << str.atoms.at(k).fpos(1) << " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
        k++;
      }
    }
  }
  if (str.num_each_type.size() == 2) {
    for (size_t i = 0; i < str.num_each_type.size(); i++) {
      for (int j = 1; j <= str.num_each_type[i]; j++) {
        if (str.num_each_type.at(0) > str.num_each_type.at(1)) {
          cout << el_names[1 - i] << j << " " << str.atoms.at(k).fpos(1) << " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
        } else {
          cout << el_names[i] << j << " " << str.atoms.at(k).fpos(1) << " " << str.atoms.at(k).fpos(2) << " " << str.atoms.at(k).fpos(3) << endl;
        }
        k++;
      }
    }
  }
  // oss <<  num_each_type.at(i) << " ";
  //  int str.atoms.size()=cpos.size();
  // for(i=0;i<str.atoms.size();i++)
  //  if(str.coord_flag==0)
  // oss << " " << str.atoms.at(i).fpos(1)<< " " << str.atoms.at(i).fpos(2) << " " << stry.atoms.at(i).fpos(3) << endl;
  //  printf("CALC ADDSYM ");
  cout << "CALC ADDSYM "; // << endl;
  if (P_EQUAL) {
    cout << "EQUAL ";
  }
  if (P_EXACT) {
    cout << "EXACT ";
  }
  if (!aurostd::isequal(P_ang, PLATON_TOLERANCE_ANGLE, 1.0e-6) || !aurostd::isequal(P_d1, PLATON_TOLERANCE_D1, 1.0e-6) || !aurostd::isequal(P_d1, PLATON_TOLERANCE_D2, 1.0e-6) ||
      !aurostd::isequal(P_d1, PLATON_TOLERANCE_D3, 1.0e-6)) {
    //    printf("%14.10f %14.10f %14.10f %14.10f ",P_ang,P_d1,P_d2,P_d3);  printf("\n");
    cout << P_ang << " " << P_d1 << " " << P_d2 << " " << P_d3 << " "; //
  }
  cout << endl;
  return cout.str();
}

// ***************************************************************************
// Function DecorateWithElements()
// ***************************************************************************
void xstructure::DecorateWithElements() {
  // Apply an element to each atom type.
  // Elements are first alphabetized to follow the AFLOW convention

  // elements need to be alphabetic for AFLOW
  deque<string> elements;
  for (size_t i = 0; i < velement.size(); i++) {
    elements.push_back(velement[i].symbol);
  } // from xelement
  std::stable_sort(elements.begin(), elements.end());

  // update species and atom names;
  SetSpecies(elements);
}

// ***************************************************************************
// xstructure::DecorateWithFakeElements() //DX20200724
// ***************************************************************************
void xstructure::DecorateWithFakeElements() {
  // Apply a fake letter to each atom type.
  // Using letters (not elements) to avoid confusion with real materials
  // i.e., prototype vs material
  // In the case of compounds with more
  // than 26 species it is necessary to add more characters to this string

  // get fake elements
  const vector<string> fake_elements = pflow::getFakeElements(num_each_type.size());

  // update species atom names;
  SetSpecies(aurostd::vector2deque(fake_elements));
}

// ***************************************************************************
// Function SpaceGroup
// ***************************************************************************
string xstructure::platon2sg(bool P_EQUAL, bool P_EXACT, double P_ang, double P_d1, double P_d2, double P_d3) {
  xstructure str = *this;
  stringstream aus;
  // string directory="/tmp/_aflow_platon_"+XHOST.ostrPID.str()+"_"+XHOST.ostrTID.str();  // dont change what works //CO20200502 - threadID
  const string directory = XHOST.tmpfs + "/_aflow_platon_" + XHOST.ostrPID.str() + "_" + XHOST.ostrTID.str();
  // dont change what works  //CO20200502 - threadID
  const string file = directory + "/aflow_platon_";
  const string file_spf = file + ".spf";
  const string file_out = file + ".out";
  vector<string> space_group;
  aurostd::DirectoryMake(directory);
  str.DecorateWithElements(); // DX20200727 - FakeNames() -> DecorateWithElements();
  aus << str.platon2print(P_EQUAL, P_EXACT, P_ang, P_d1, P_d2, P_d3);
  aurostd::stringstream2file(aus, file_spf);
  aus.clear();
  aus.str(std::string());
  aus << "cd " << directory << endl;
  aus << "platon -o -c " << file_spf << R"( | grep "Space Group" | head -1 | sed "s/Space Group //g"  | sed "s/No:/#/g" | sed "s/,/#/g" | sed "s/ //g" > )" << file_out << endl;
  aurostd::execute(aus);
  aurostd::string2tokens(aurostd::file2string(file_out), space_group, "#");
  aurostd::RemoveDirectory(directory);
  spacegroup = space_group[0] + " #" + space_group[1];
  spacegrouplabel = space_group[0];
  spacegroupnumber = aurostd::string2utype<int>(space_group[1]);
  spacegroupoption = "FIX: to extract from platon";
  is_spacegroup_platon = true;
  return spacegroup;
}

string xstructure::findsym2sg(double tolerance) {
  xstructure str = *this;
  // Read in input file.
  stringstream cout;
  const string findsym_in = aurostd::TmpFileCreate("findsym.in");
  cout << str.findsym2print(tolerance);
  aurostd::stringstream2file(cout, findsym_in);
  vector<string> vlines;
  vector<string> tokens;
  FINDSYM::Write("data_space.txt", "./");
  FINDSYM::Write("data_wyckoff.txt", "./");
  aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("findsym") + " < " + findsym_in), vlines);
  FROZSL::Delete("data_space.txt", "./");
  FROZSL::Delete("data_wyckoff.txt", "./");
  aurostd::RemoveFile(findsym_in);
  aurostd::RemoveFile("./findsym.log");
  for (size_t i = 0; i < vlines.size(); i++) {
    if (aurostd::substring2bool(vlines[i], "Space Group")) {
      aurostd::string2tokens(vlines[i], tokens);
    }
  }
  //    for(size_t i=0;i<tokens.size();i++) cerr << tokens.at(i) << endl;
  spacegroup = tokens.at(4) + " #" + tokens.at(2);
  spacegrouplabel = tokens.at(4);
  spacegroupnumber = aurostd::string2utype<int>(tokens.at(2));
  spacegroupoption = "FIX: to extract from findsym";
  is_spacegroup_findsym = true;
  return spacegroup;
}

// ***************************************************************************
// xstructure::findsym2execute
// ***************************************************************************
// This funtion executes findsys
string xstructure::findsym2execute(double tolerance) {
  xstructure str = *this;
  // Read in input file.
  stringstream cout;
  const string findsym_in = aurostd::TmpFileCreate("findsym.in");
  cout << str.findsym2print(tolerance);
  aurostd::stringstream2file(cout, findsym_in);
  const vector<string> vlines;
  const vector<string> tokens;
  FINDSYM::Write("data_space.txt", "./");
  FINDSYM::Write("data_wyckoff.txt", "./");
  string out = aurostd::execute2string(XHOST.command("findsym") + " < " + findsym_in);
  FROZSL::Delete("data_space.txt", "./");
  FROZSL::Delete("data_wyckoff.txt", "./");
  aurostd::RemoveFile(findsym_in);
  aurostd::RemoveFile("./findsym.log");
  return out;
}

// ***************************************************************************
// xstructure::findsym2print
// ***************************************************************************
// This funtion prints out structural data in a format for findsys
string xstructure::findsym2print(double tolerance) {
  xstructure sstr = *this;
  stringstream cout;
  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr.ReScale(1.0);
  // Print out data in Findsym format
  cout << sstr.title << endl;
  cout << tolerance << "                               accuracy (0=highest=1e-6) " << endl;
  cout << "2 Form of lattice parameters: to be entered as lengths and angles " << endl;
  cout << sstr.a << " " << sstr.b << " " << sstr.c << " " << sstr.alpha << " " << sstr.beta << " " << sstr.gamma << " a, b, c, alpha, beta, gamma " << endl;
  cout << "1 Vectors defining the unit cell" << endl;
  cout << "1.0 0.0 0.0                   unit cell in function of lattice parameter" << endl;
  cout << "0.0 1.0 0.0                   unit cell in function of lattice parameter" << endl;
  cout << "0.0 0.0 1.0                   unit cell in function of lattice parameter" << endl;
  cout << sstr.atoms.size() << "                            number of atoms in the primitive unit cell" << endl;
  for (size_t i = 0; i < sstr.atoms.size(); i++) {
    cout << sstr.atoms[i].type + 1 << " ";
  }
  cout << endl;
  for (size_t i = 0; i < sstr.atoms.size(); i++) {
    //  if(coord_flag==0)
    cout << " " << sstr.atoms[i].fpos(1) << " " << sstr.atoms[i].fpos(2) << " " << sstr.atoms[i].fpos(3) << endl;
    //  if(coord_flag==1) oss << " " << sstr.atoms[i].cpos(1) << " " << sstr.atoms[i].cpos(2) << " " << sstr.atoms[i].cpos(3) << endl;
  }
  //  oss << endl;
  return cout.str();
}

// Stefano Curtarolo MIT 2002 - DUKE 2013

// ***************************************************************************
// Rotate
// ***************************************************************************
// Rotation is done around the origin of the structure.
// Define origin q, a point p, and a rotation of p around
// q (R_q(p)).  Then R_q(p)=R_0(p-q)+q=R_0(p)-R_0(q)+q.
// We can evalutate the R_0 terms by simply mulitplying by
// a matrix.  Then we must add q to all the final cartesian
// positions.

// ---------------------------------------------------------------------------
// returns new xstructure (makes a copy)
xstructure Rotate(const xstructure& a, const xmatrix<double>& rm) {
  xstructure xstr_rotated = a;
  xstr_rotated.Rotate(rm);
  return xstr_rotated;
}

// ---------------------------------------------------------------------------
// modifies in-place (efficient)
void xstructure::Rotate(const xmatrix<double>& rm) {
  const bool LDEBUG = (false || XHOST.DEBUG); // CO20190520
  if (LDEBUG) { // CO20190520
    cerr << __AFLOW_FUNC__ << " (*this)=" << endl;
    cerr << (*this) << endl; // CO20190520
    cerr << __AFLOW_FUNC__ << " (*this).origin=" << (*this).origin << endl; // CO20190520
    cerr << __AFLOW_FUNC__ << " rm=" << endl;
    cerr << rm << endl; // ME20200204
  }
  if (aurostd::isidentity(rm)) {
    return; // ME20200204 - no need to go through all the motions for identity matrix
  }
  // Get R_0(p) for all cartesian positions.
  xmatrix<double> nlattice(3, 3);
  nlattice = trasp((*this).lattice);
  (*this).lattice = trasp(rm * nlattice);
  (*this).FixLattices(); // CO20190409 - so we don't need to keep redefining f2c/c2f
  const xmatrix<double>& f2c = (*this).f2c; // CO20190520
  const xmatrix<double>& c2f = (*this).c2f; // CO20190520
  const size_t natoms = (*this).atoms.size(); // DX+ME20210111 - set variable to optimize for-loops
  // Get R_0(q)
  xvector<double> r_orig(3);
  r_orig = rm * (*this).origin;
  // Assign new cartesian positions
  // Get all the direct coords.
  for (uint ia = 0; ia < natoms; ia++) { // DX+ME20210111 consolidate into one for-loop
    (*this).atoms[ia].cpos = f2c * (*this).atoms[ia].fpos;
    (*this).atoms[ia].cpos += (-r_orig + (*this).origin);
    (*this).atoms[ia].fpos = c2f * (*this).atoms[ia].cpos;
  }
  return;
}

// ***************************************************************************
// GetRotationMatrix // moved from pflow
// ***************************************************************************
// This gets a rotation matrix from 3 angles assumed
// to represent a rotation around x, then y, then z.
// Angles are assumed to be in radians.
aurostd::matrix<double> GetRotationMatrix(const vector<double>& angles) {
  // CO20200404 pflow::matrix()->aurostd::matrix()
  //  Sin and cos.
  vector<double> sn(3, 0.0);
  vector<double> cs(3, 0.0);
  for (int ic = 0; ic < 3; ic++) {
    sn[ic] = sin(angles[ic]);
    cs[ic] = cos(angles[ic]);
  }
  // Set rotation matrix (do x, then y, then z rotation).
  aurostd::matrix<double> xm(3, 3);
  pflow::VVset(xm, 0.0); // CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> ym(3, 3);
  pflow::VVset(ym, 0.0); // CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> zm(3, 3);
  pflow::VVset(zm, 0.0); // CO20200404 pflow::matrix()->aurostd::matrix()

  xm[0][0] = 1;
  xm[1][1] = cs[0];
  xm[1][2] = -sn[0];
  xm[2][1] = sn[0];
  xm[2][2] = cs[0];

  ym[0][0] = cs[1];
  ym[0][2] = sn[1];
  ym[1][1] = 1;
  ym[2][0] = -sn[1];
  ym[2][2] = cs[1];

  zm[0][0] = cs[2];
  zm[0][1] = -sn[2];
  zm[1][0] = sn[2];
  zm[1][1] = cs[2];
  zm[2][2] = 1;

  aurostd::matrix<double> rm; // CO20200404 pflow::matrix()->aurostd::matrix()
  rm = pflow::MMmult(zm, pflow::MMmult(ym, xm));
  return rm;
}

// ***************************************************************************
// RotateStrVec // moved from pflow
// ***************************************************************************
// This rotates each structure in the vstr.
// The rotation goes from an initial to a final set of angles
// in steps of the primitive rotation.
// The primitive rotation is the rotation around x, y, then z
// by the amount of change given in rot divided by the number
// of structures - 1.  All frames are initially rotated
// according to the initial rotation angles.  The rotation is
// then done as one primitive rotation per frame.
void RotateStrVec(vector<xstructure>& vstr, const vector<double>& rot) {
  // Get initial rotation matrix.
  vector<double> angles(3);
  for (int ic = 0; ic < 3; ic++) {
    angles[ic] = rot[2 * ic];
    angles[ic] = angles[ic] * TWOPI / 360.0;
  }
  const aurostd::matrix<double> irm = GetRotationMatrix(angles); // CO20200404 pflow::matrix()->aurostd::matrix()

  // get primitive rotation matrix.
  int s = vstr.size() - 1;
  if (s < 1) {
    s = 1;
  }
  for (int ic = 0; ic < 3; ic++) {
    angles[ic] = (rot[2 * ic + 1] - rot[2 * ic]) / (s);
    angles[ic] = angles[ic] * TWOPI / 360.0;
  }
  const aurostd::matrix<double> prm = GetRotationMatrix(angles); // CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> rm = irm; // CO20200404 pflow::matrix()->aurostd::matrix()
  xmatrix<double> xprm(3, 3);
  xprm = aurostd::matrix2xmatrix(prm); // CO20200404 pflow::matrix()->aurostd::matrix()
  xmatrix<double> xrm(3, 3);
  xrm = aurostd::matrix2xmatrix(rm); // CO20200404 pflow::matrix()->aurostd::matrix()
  for (int is = 0; is < (int) vstr.size(); is++) {
    //    xrm=aurostd::matrix2xmatrix(rm); //CO20200404 pflow::matrix()->aurostd::matrix()
    //   vstr[is]=Rotate(vstr[is],xrm);
    vstr[is].Rotate(aurostd::matrix2xmatrix(rm));
    // CO20200404 pflow::matrix()->aurostd::matrix() //DX20210127 - do not make copy of xstructure
    rm = pflow::MMmult(prm, rm);
  }
}

// ***************************************************************************
//  Function GetLTCell
// ***************************************************************************
// This function operates a linear transform on a cell.
xstructure GetLTCell(const xmatrix<double>& lt, const xstructure& str) {
  xstructure nstr(str);
  nstr.lattice = str.lattice * trasp(lt);
  for (int iat = 0; iat < (int) nstr.atoms.size(); iat++) {
    //    nstr.atoms[iat].cpos=str.atoms[iat].cpos*trasp(lt);
    //  nstr.atoms[iat].fpos=C2F(nstr.lattice,nstr.atoms[iat].cpos);
    nstr.atoms[iat] = F2C(nstr.lattice, nstr.atoms[iat]);
  }
  return nstr;
}

// ***************************************************************************
//  Function GetLTFVCell
// ***************************************************************************
// This function operates a linear transform on a cell.  Rotates all vectors
// by phi around nvec.  Formula from Goldstein, 1980 (Actually used version on
// Mathematic web site).

xstructure GetLTFVCell(const xvector<double>& nvec, const double phi, const xstructure& str) {
  //  xmatrix<double> fpos=str.fpos;
  const xmatrix<double> nlat(3, 3);
  xvector<double> nhat(3);
  xvector<double> v1(3);
  xvector<double> v2(3);
  xvector<double> v3(3);
  xvector<double> r(3);
  nhat = nvec / aurostd::modulus(nvec);
  double tmp;
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  for (int i = 1; i <= 3; i++) {
    r = str.lattice(i);
    const double nhatpr = scalar_product(nhat, r);
    v1 = cphi * r;
    tmp = nhatpr * (1.0 - cphi);
    v2 = tmp * nhat;
    v3 = vector_product(r, nhat);
    v3 = sphi * v3;
    for (int j = 1; j <= 3; j++) {
      nlat(i, j) = nlat(i, j) + v1(j);
    }
    for (int j = 1; j <= 3; j++) {
      nlat(i, j) = nlat(i, j) + v2(j);
    }
    for (int j = 1; j <= 3; j++) {
      nlat(i, j) = nlat(i, j) + v3(j);
    }
  }
  xstructure nstr = str;
  nstr.lattice = nlat;
  for (int iat = 0; iat < (int) nstr.atoms.size(); iat++) {
    nstr.atoms[iat].fpos = str.atoms.at(iat).fpos; // *trasp(lt);
    nstr.atoms[iat].cpos = F2C(nstr.lattice, nstr.atoms[iat].fpos);
  }
  return nstr;
}

// ***************************************************************************
//  Function ShiftPos ShiftCpos ShiftFpos
// ***************************************************************************
// This function shifts the position of the atoms !
// DX20201215 - added in-place methods

// make a copy
xstructure ShiftPos(const xstructure& a, const xvector<double>& shift, bool is_frac) {
  // DX20210113 - "int flag" to "bool is_frac"
  xstructure b(a);
  b.ShiftPos(shift, is_frac);
  return b;
}

// modify in-place //DX2021215
void xstructure::ShiftPos(const xvector<double>& shift, bool is_frac) {
  // Cartesian shift
  if (is_frac == false) {
    (*this).ShiftCPos(shift);
  } // DX20210113 - use function, reduce code
  // Direct coords shift
  else {
    (*this).ShiftFPos(shift);
  } // DX20210113 - use function, reduce code
}

// make a copy
xstructure ShiftCPos(const xstructure& a, const xvector<double>& shift) {
  xstructure b(a);
  b.ShiftCPos(shift);
  return b;
}

// modify in-place //DX2021215
void xstructure::ShiftCPos(const xvector<double>& shift) {
  const size_t natoms = (*this).atoms.size(); // DX20210111 - make variable to optimize for-loops
  const xmatrix<double> c2f = (*this).scale * inverse(trasp((*this).lattice)); // DX+ME20210111
  (*this).coord_flag = true;
  for (uint ia = 0; ia < natoms; ia++) {
    (*this).atoms[ia].cpos = (*this).atoms[ia].cpos + shift;
    (*this).atoms[ia].fpos = c2f * (*this).atoms[ia].cpos; // DX+ME20210111 - use pre-calculated c2f, optimize
  }
}

// make a copy
xstructure ShiftFPos(const xstructure& a, const xvector<double>& shift) {
  xstructure b(a);
  b.ShiftCPos(shift);
  return b;
}

// modify in-place //DX2021215
void xstructure::ShiftFPos(const xvector<double>& shift) {
  const size_t natoms = (*this).atoms.size(); // DX20210111 - make variable to optimize for-loops
  const xmatrix<double> f2c = (*this).scale * trasp((*this).lattice); // DX+ME20210111
  (*this).coord_flag = false;
  for (uint ia = 0; ia < natoms; ia++) {
    (*this).atoms[ia].fpos = (*this).atoms[ia].fpos + shift;
    (*this).atoms[ia].cpos = f2c * (*this).atoms[ia].fpos; // DX+ME20210111 - use pre-calculated f2c, optimize
  }
}

// **************************************************************************
// MaxStructureLattice and MinStructureLattice
// **************************************************************************
double MaxStructureLattice(const xstructure& str) {
  return max(aurostd::modulus(str.lattice(1)), aurostd::modulus(str.lattice(2)), aurostd::modulus(str.lattice(3)));
}

double MinStructureLattice(const xstructure& str) {
  return min(aurostd::modulus(str.lattice(1)), aurostd::modulus(str.lattice(2)), aurostd::modulus(str.lattice(3)));
}

// **************************************************************************
// AtomCDisp
// **************************************************************************
// This function returns the Cartesian displacement vector at2-at1 between two atoms.
xvector<double> AtomCDisp(const _atom& at1, const _atom& at2) {
  xvector<double> diff(3);
  diff = at2.cpos - at1.cpos;
  return diff;
}

// **************************************************************************
// AtomDist
// **************************************************************************
// This function returns the distance between two atoms.
double AtomDist(const xstructure& str, const _atom& atom1, const _atom& atom2) { // with structure
  xvector<double> pos1(3);
  xvector<double> pos2(3);
  // contains shift WITH the ijk lattice and origin of structure
  pos1 = str.scale * (atom1.cpos + atom1.ijk(1) * str.lattice(1) + atom1.ijk(2) * str.lattice(2) + atom1.ijk(3) * str.lattice(3)) + str.origin;
  pos2 = str.scale * (atom2.cpos + atom2.ijk(1) * str.lattice(1) + atom2.ijk(2) * str.lattice(2) + atom2.ijk(3) * str.lattice(3)) + str.origin;
  return aurostd::modulus(pos1 - pos2);
}

double AtomDist(const _atom& at1, const _atom& at2) { // without structure
  return aurostd::modulus(AtomCDisp(at1, at2));
}

// **************************************************************************
// SameAtom
// **************************************************************************
#define _atom_eps_tol_ 0.01

bool SameAtom(const xstructure& str, const _atom& atom1, const _atom& atom2) {
  if (AtomDist(str, atom1, atom2) < _atom_eps_tol_ && atom1.type == atom2.type) {
    return true;
  }
  return false;
}

bool SameAtom(const _atom& atom1, const _atom& atom2) {
  if (AtomDist(atom1, atom2) < _atom_eps_tol_ && atom1.type == atom2.type) {
    return true;
  }
  return false;
}

// **************************************************************************
// DifferentAtom
// **************************************************************************
bool DifferentAtom(const xstructure& str, const _atom& atom1, const _atom& atom2) {
  if (AtomDist(str, atom1, atom2) > _atom_eps_tol_ || atom1.type != atom2.type) {
    return true;
  }
  return false;
}

// **************************************************************************
// GetDistMatrix
// **************************************************************************
// CO20171024
xmatrix<double> GetDistMatrix(const xstructure& aa) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  xstructure a(aa);
  a.ReScale(1.0);
  xstructure xstr_cluster = a;
  // int dim=aurostd::max(LatticeDimensionSphere(a.lattice,RadiusSphereLattice(a.lattice))); //OVERKILL
  GenerateGridAtoms(xstr_cluster, RadiusSphereLattice(a.lattice)); // dim,dim,dim);            //OVERKILL
  vector<vector<int>> atom_types;
  atom_types.emplace_back(0);
  atom_types.back().push_back(0);
  bool found;
  uint indx_found;
  for (size_t i = 1; i < a.atoms.size(); i++) {
    found = false;
    for (size_t j = 0; j < atom_types.size() && !found; j++) {
      if (atom_types[j].empty()) {
        cerr << "GetDistMatrix: atom_types indices are compromised!" << endl;
        const xmatrix<double> dummy; // dummy
        return dummy;
      } // should never come up
      if (a.atoms[i].type == a.atoms[atom_types[j][0]].type) {
        found = true;
        indx_found = j;
      }
    }
    if (!found) {
      atom_types.emplace_back(0);
      indx_found = atom_types.size() - 1;
    }
    atom_types[indx_found].push_back(i);
  }
  if (!a.num_each_type.empty()) {
    if (atom_types.size() != a.num_each_type.size()) {
      cerr << "GetDistMatrix: odd count of atom types, it's probably poorly sorted" << endl;
      const xmatrix<double> dummy; // dummy
      return dummy;
    }
  }
  if (LDEBUG) {
    for (size_t it1 = 0; it1 < atom_types.size(); it1++) {
      cerr << "GetDistMatrix: atom_types[" << it1 << "][0]=" << atom_types[it1][0] << endl;
    }
  }
  const xmatrix<double> distsij(1, 1, atom_types.size(), atom_types.size());
  for (size_t it1 = 0; it1 < atom_types.size(); it1++) {
    for (size_t it2 = 0; it2 < atom_types.size(); it2++) {
      distsij(it1 + 1, it2 + 1) = distsij(it2 + 1, it1 + 1) = AUROSTD_MAX_DOUBLE;
    }
  }
  uint atom1 = 0;
  uint atom2 = 0;
  uint it1 = 0;
  uint it2 = 0;
  uint ia1 = 0;
  uint ia2 = 0;
  double distij;
  double min_dist;
  for (it1 = 0; it1 < atom_types.size(); it1++) { // type 1 (must be in this order)
    for (it2 = it1; it2 < atom_types.size(); it2++) { // type 2 (must be in this order)
      if (LDEBUG) {
        cerr << "GetDistMatrix: finding min dist between itype=" << it1 << " and itype=" << it2 << endl;
      }
      min_dist = AUROSTD_MAX_DOUBLE; // reset min_dist
      for (ia1 = 0; ia1 < atom_types[it1].size(); ia1++) { // must go through all atoms of same type
        atom1 = xstr_cluster.grid_atoms_pc2scMap[atom_types[it1][ia1]]; // get respective index in cluster
        for (ia2 = 0; ia2 < atom_types[it2].size(); ia2++) { // must go through all atoms of the same type
          for (atom2 = 0; atom2 < (uint) xstr_cluster.grid_atoms_number; atom2++) { // go through all atoms of the cluster
            if (atom1 != atom2 && a.atoms[atom_types[it2][ia2]].type == xstr_cluster.grid_atoms[atom2].type) {
              // cannot be same index (dist=0), and must be the types we want
              distij = AtomDist(xstr_cluster.grid_atoms[atom1], xstr_cluster.grid_atoms[atom2]); // distance
              if (false && LDEBUG) {
                cerr << "GetDistMatrix: trying dist=" << distij << " between atom[" << atom1 << ",type=" << xstr_cluster.grid_atoms[atom1].type << "] and atom[" << atom2
                     << ",type=" << xstr_cluster.grid_atoms[atom2].type << "]" << endl;
              }
              if (distij < min_dist) { // grab min_dist
                distsij(it1 + 1, it2 + 1) = distsij(it2 + 1, it1 + 1) = distij;
                min_dist = distij;
                if (LDEBUG) {
                  cerr << "GetDistMatrix: found new min dist=" << min_dist << " between atom[" << atom1 << ",type=" << xstr_cluster.grid_atoms[atom1].type << "] and atom[" << atom2
                       << ",type=" << xstr_cluster.grid_atoms[atom2].type << "]" << endl;
                }
              }
            }
          }
        }
      }
    }
  }
  return distsij;
}

// **************************************************************************
// GetNBONDXX
// **************************************************************************
// CO20171024
vector<double> GetNBONDXX(const xstructure& a) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  const xmatrix<double> distsij = GetDistMatrix(a);
  vector<double> dists;
  for (uint it1 = 1; it1 < (uint) distsij.rows + 1; it1++) {
    for (uint it2 = it1; it2 < (uint) distsij.cols + 1; it2++) {
      dists.push_back(distsij(it1, it2));
    }
  }
  if (LDEBUG) {
    cerr << "GetNBONDXX: DIST_MATRIX" << endl;
    cerr << distsij << endl;
    cerr << "GetNBONDXX: MIN_DIST=" << SYM::minimumDistance(a) << endl;
  }
  return dists;
}

// **************************************************************************
// GenerateGridAtoms
// **************************************************************************
// make grid of atoms !
int xstructure::GenerateGridAtoms(double radius) {
  return GenerateGridAtoms(LatticeDimensionSphere((*this), radius));
}
// radius is not normalized over the scale
int xstructure::GenerateGridAtoms(int d) {
  return GenerateGridAtoms(-d, d, -d, d, -d, d);
}
int xstructure::GenerateGridAtoms(int d1, int d2, int d3) {
  return GenerateGridAtoms(-d1, d1, -d2, d2, -d3, d3);
}

int xstructure::GenerateGridAtoms(const xvector<int>& dims) {
  return GenerateGridAtoms(-dims(1), dims(1), -dims(2), dims(2), -dims(3), dims(3));
} // CO20200912
int xstructure::GenerateGridAtoms(int i1, int i2, int j1, int j2, int k1, int k2) { // DX20191218 - [NEW]
  const bool LDEBUG = (false || XHOST.DEBUG); // CO20190520
  if (LDEBUG) { // CO20190520
    cerr << __AFLOW_FUNC__ << " str=" << endl;
    cerr << (*this) << endl; // CO20190520
    cerr << __AFLOW_FUNC__ << " i=" << i1 << ":" << i2 << endl; // CO20190520
    cerr << __AFLOW_FUNC__ << " j=" << j1 << ":" << j2 << endl; // CO20190520
    cerr << __AFLOW_FUNC__ << " k=" << k1 << ":" << k2 << endl; // CO20190520
  } // CO20190520
  // same scale as before
  grid_atoms.clear();
  grid_atoms_sc2pcMap.clear();
  grid_atoms_pc2scMap.clear();
  _atom atom;
  BringInCell(); // are INCELL.
  // xvector<double> a1(3),a2(3),a3(3);                     // a1,a2,a3 are the rows of the lattice matrix
  // a1=lattice(1);a2=lattice(2);a3=lattice(3); // a1,a2,a3 are the rows of the lattice matrix
  const xvector<double>& a1 = lattice(1); // CO20190520 - no need to make copies
  const xvector<double>& a2 = lattice(2); // CO20190520 - no need to make copies
  const xvector<double>& a3 = lattice(3); // CO20190520 - no need to make copies
  // DX20190709 - calculate and store once = speed - START
  vector<xvector<double>> l1;
  vector<xvector<double>> l2;
  vector<xvector<double>> l3;
  vector<int> a_index;
  vector<int> b_index;
  vector<int> c_index;
  for (int i = i1; i <= i2; i++) {
    l1.push_back(i * a1);
    a_index.push_back(i);
  }
  for (int j = j1; j <= j2; j++) {
    l2.push_back(j * a2);
    b_index.push_back(j);
  }
  for (int k = k1; k <= k2; k++) {
    l3.push_back(k * a3);
    c_index.push_back(k);
  }
  // DX20191218 - calculate and store once = speed - END

  // resize vectors - DX20191122
  const size_t num_grid_atoms = atoms.size() * l1.size() * l2.size() * l3.size();
  grid_atoms.resize(num_grid_atoms);
  grid_atoms_pc2scMap.resize(atoms.size()); // DX20191218 - should be the size of the the primitive cell, not grid
  grid_atoms_sc2pcMap.resize(num_grid_atoms);

  uint grid_atom_count = 0; // keep track of index - DX20191122

  for (size_t iat = 0; iat < atoms.size(); iat++) {
    // grid_atoms.push_back(atoms[iat]);  // put first the unit cell ! //DX20190709 - at to [] = speed increase
    // grid_atoms_pc2scMap.push_back(grid_atoms.size()-1); //CO20171025
    // grid_atoms_sc2pcMap.push_back(iat); //CO20171025
    grid_atoms[grid_atom_count] = atoms[iat]; // put first the unit cell ! //DX20190709 - at to [] = speed increase
    grid_atoms_pc2scMap[grid_atom_count] = iat;
    // DX20191218 - use the index of the primitive cell not the running count of grid_atoms since it was resized
    grid_atoms_sc2pcMap[grid_atom_count] = iat; // CO20171025
    grid_atom_count++; // DX20191122
  }

  xvector<double> a_component(3);
  xvector<double> ab_component(3);
  xvector<double> abc_component(3);
  // DX+ME20191107 - define outside loop (speed increase)
  const size_t natoms = atoms.size(); // DX20191107 - initialize natoms outside loop (speed increase)
  for (size_t i = 0; i < l1.size(); i++) {
    a_component = l1[i]; // DX : i*lattice(1)
    for (size_t j = 0; j < l2.size(); j++) {
      ab_component = a_component + l2[j]; // DX : i*lattice(1) + j*lattice(2)
      for (size_t k = 0; k < l3.size(); k++) {
        // DX20191218 [WRONG INDICES] if(i!=0 || j!=0 || k!=0)
        if (a_index[i] != 0 || b_index[j] != 0 || c_index[k] != 0) // DX20191218
        { // CO20200106 - patching for auto-indenting
          abc_component = ab_component + l3[k]; // DX : i*lattice(1) + j*lattice(2) + k*lattice(3)
          for (size_t iat = 0; iat < natoms; iat++) { // DX20191107 - replace atoms.size() with natoms
            atom = atoms[iat]; // DX20190709 - at to [] = speed increase
            atom.isincell = false; // these are OUT OF CELL
            atom.cpos += abc_component; // DX20190709 - at to [] = speed increase //CO20191127
            atom.fpos[1] += a_index[i]; // DX20190709 - at to [] = speed increase //CO20191127
            atom.fpos[2] += b_index[j]; // DX20190709 - at to [] = speed increase //CO20191127
            atom.fpos[3] += c_index[k]; // DX20190709 - at to [] = speed increase //CO20191127
            grid_atoms[grid_atom_count] = atom;
            grid_atoms_sc2pcMap[grid_atom_count] = iat; // CO20171025
            grid_atom_count++; // DX20191122
          }
        }
      }
    }
  }
  // DX20200320 - moved outside of loop so that LDEBUG boolean is not checked every time (speed increase when grid atoms is large)
  if (LDEBUG) { // CO20190520
    for (size_t i = 0; i < grid_atoms.size(); i++) {
      cerr << __AFLOW_FUNC__ << " grid_atoms[" << i << "].cpos=" << grid_atoms[i].cpos << endl; // CO20190520
      cerr << __AFLOW_FUNC__ << " grid_atoms[" << i << "].fpos=" << grid_atoms[i].fpos << endl; // CO20190520
      cerr << __AFLOW_FUNC__ << " grid_atoms[" << i << "]=" << grid_atoms[i] << endl; // DX20191218
      cerr << __AFLOW_FUNC__ << " grid_atoms_sc2pcMap[" << i << "]=" << grid_atoms_sc2pcMap[i] << endl; // DX20191218
    } // CO20190520
  }
  grid_atoms_calculated = true;
  grid_atoms_dimsL[1] = i1;
  grid_atoms_dimsL[2] = j1;
  grid_atoms_dimsL[3] = k1;
  grid_atoms_dimsH[1] = i2;
  grid_atoms_dimsH[2] = j2;
  grid_atoms_dimsH[3] = k2;
  grid_atoms_number = grid_atoms.size();
  //  for(size_t i=0;i<grid_atoms.size();i++)
  //   cerr << grid_atoms.at(i).cpos << endl;
  // cerr << grid_atoms.size() << endl;
  return grid_atoms.size();
}

int GenerateGridAtoms(xstructure& str) {
  return str.GenerateGridAtoms(-1, 1, -1, 1, -1, 1);
}
int GenerateGridAtoms(xstructure& str, double radius) {
  return str.GenerateGridAtoms(radius);
} // CO20200912 - double
int GenerateGridAtoms(xstructure& str, int d) {
  return str.GenerateGridAtoms(d);
}
int GenerateGridAtoms(xstructure& str, int d1, int d2, int d3) {
  return str.GenerateGridAtoms(d1, d2, d3);
}

int GenerateGridAtoms(xstructure& str, int i1, int i2, int j1, int j2, int k1, int k2) {
  return str.GenerateGridAtoms(i1, i2, j1, j2, k1, k2);
} // DX20191218
int GenerateGridAtoms(xstructure& str, const xvector<int>& dims) {
  return str.GenerateGridAtoms(dims);
} // CO20200912

// **************************************************************************
// GenerateLIJK table stuff
// **************************************************************************
namespace aurostd { // INT
  class _ssort_int_value0123 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const {
      if (v1[0] < v2[0]) {
        return true;
      }
      if (v1[0] == v2[0]) {
        if (v1[1] > v2[1]) {
          return true;
        }
        if (v1[1] == v2[1]) {
          if (v1[2] > v2[2]) {
            return true;
          }
          if (v1[2] == v2[2]) {
            if (v1[3] > v2[3]) {
              return true;
            }
          }
        }
      }
      return false;
    }
  };
} // namespace aurostd

int xstructure::GenerateLIJK(double radius) {
  lijk_table.clear();
  xvector<double> a1(3);
  xvector<double> a2(3);
  xvector<double> a3(3); //    a1,a2,a3 are the rows of the lattice matrix
  a1 = lattice(1);
  a2 = lattice(2);
  a3 = lattice(3); // a1,a2,a3 are the rows of the lattice matrix
  vector<int> _lijk(4);
  const xvector<int> int_ijk(3);
  xvector<double> cpos(3);
  const xvector<double> fpos(3);
  xvector<int> dims(3);
  dims = LatticeDimensionSphere(lattice, radius); // for(int i=0;i<=3;i++) dims[i]+=1;
  lijk_dims[1] = dims[1];
  lijk_dims[2] = dims[2];
  lijk_dims[3] = dims[3];
  vector<vector<int>> _temp_lijk_table;
  _temp_lijk_table.clear();

  for (int i = -lijk_dims[1]; i <= lijk_dims[1]; i++) {
    for (int j = -lijk_dims[2]; j <= lijk_dims[2]; j++) {
      for (int k = -lijk_dims[3]; k <= lijk_dims[3]; k++) {
        cpos = ((double) i) * a1 + ((double) j) * a2 + ((double) k) * a3;
        // 	if(aurostd::modulus(cpos)<=radius)  // I put them all since I will not be scanning through.
        {
          _lijk[0] = (int) (10000000.0 * sqrt(aurostd::modulus(cpos))); // for the next sorting
          _lijk[1] = i;
          _lijk[2] = j;
          _lijk[3] = k;
          _temp_lijk_table.push_back(_lijk);
        }
      }
    }
  }
  sort(_temp_lijk_table.begin(), _temp_lijk_table.end(), aurostd::_ssort_int_value0123());

  for (size_t i = 0; i < _temp_lijk_table.size(); i++) {
    _temp_lijk_table[i][0] = i; // relabel
    int_ijk[1] = _temp_lijk_table[i][1];
    int_ijk[2] = _temp_lijk_table[i][2];
    int_ijk[3] = _temp_lijk_table[i][3];
    lijk_table.push_back(int_ijk);
  }

  for (size_t i = 0; i < lijk_table.size(); i++) {
    cpos = ((double) lijk_table[i][1]) * a1 + ((double) lijk_table[i][2]) * a2 + ((double) lijk_table[i][3]) * a3;
    fpos[1] = lijk_table[i][1];
    fpos[2] = lijk_table[i][2];
    fpos[3] = lijk_table[i][3];
    lijk_cpos.push_back(cpos);
    lijk_fpos.push_back(fpos);
  }

  lijk_calculated = true;
  return lijk_table.size();
}

// **************************************************************************
// LIJK table: L => IJK
// **************************************************************************
void l2ijk(const xstructure& str, const int& l, int& i, int& j, int& k) {
  if (l < 0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "l2ijk():", "l2ijk error: l<0", _VALUE_RANGE_);
  }
  if (l > (int) str.lijk_table.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "l2ijk():", "l>str.lijk_table.size()", _VALUE_RANGE_);
  }
  i = str.lijk_table[l][1];
  j = str.lijk_table[l][2];
  k = str.lijk_table[l][3];
}

void l2ijk(const xstructure& str, const int& l, xvector<int>& ijk) {
  if (l < 0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "l2ijk():", "l2ijk error: l<0", _VALUE_RANGE_);
  }
  if (l > (int) str.lijk_table.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "l2ijk():", "l>str.lijk_table.size()", _VALUE_RANGE_);
  }
  ijk[1] = str.lijk_table[l][1];
  ijk[2] = str.lijk_table[l][2];
  ijk[3] = str.lijk_table[l][3];
}

xvector<int> l2ijk(const xstructure& str, const int& l) {
  const xvector<int> ijk(3);
  int i;
  int j;
  int k;
  l2ijk(str, l, i, j, k);
  ijk(1) = i;
  ijk(2) = j;
  ijk(3) = k;
  return ijk;
}

// **************************************************************************
// LIJK table: IJK => L
// **************************************************************************
void ijk2l(const xstructure& str, int& l, const int& i, const int& j, const int& k) {
  l = -1;
  if (i < -str.lijk_dims(1) || i > str.lijk_dims(1)) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "ijk2l():", "i out of boundary", _VALUE_RANGE_);
  }
  if (j < -str.lijk_dims(2) || j > str.lijk_dims(2)) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "ijk2l():", "j out of boundary", _VALUE_RANGE_);
  }
  if (k < -str.lijk_dims(3) || k > str.lijk_dims(3)) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "ijk2l():", "k out of boundary", _VALUE_RANGE_);
  }
  for (size_t ll = 0; ll < str.lijk_table.size(); ll++) { // start search
    if (str.lijk_table[ll][1] == i) { // faster comparison one at a time
      if (str.lijk_table[ll][2] == j) { // faster comparison one at a time
        if (str.lijk_table[ll][3] == k) { // faster comparison one at a time
          l = ll;
        }
      }
    }
  }
  if (l < 0) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "ijk2l():", "l not found", _VALUE_RANGE_);
  }
}

void ijk2l(const xstructure& str, int& l, const xvector<int>& ijk) {
  ijk2l(str, l, ijk(1), ijk(2), ijk(3));
}

int ijk2l(const xstructure& str, const int& i, const int& j, const int& k) {
  int l;
  ijk2l(str, l, i, j, k);
  return l;
}

int ijk2l(const xstructure& str, const xvector<int>& ijk) {
  int l;
  ijk2l(str, l, ijk(1), ijk(2), ijk(3));
  return l;
}

// **************************************************************************
// input2AIMSxstr input2ABINITxstr input2QExstr input2VASPxstr
// ***************************************************************************
xstructure input2AIMSxstr(istream& input) {
  xstructure a(input, IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR)
  a.xstructure2aims();
  return a;
}

xstructure input2ABINITxstr(istream& input) {
  xstructure a(input, IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR)
  a.xstructure2abinit();
  return a;
}

xstructure input2QExstr(istream& input) {
  xstructure a(input, IOAFLOW_AUTO);
  //  if(a.iomode==IOVASP_AUTO || a.iomode==IOVASP_POSCAR || a.iomode==IOVASP_ABCCAR || a.iomode==IOVASP_WYCKCAR)
  a.xstructure2qe();
  return a;
}

xstructure input2VASPxstr(istream& input, bool vasp5) { // CO20210119 - added vasp5
  xstructure a(input, IOAFLOW_AUTO);
  //  if(a.iomode==IOQE_AUTO || a.iomode==IOQE_GEOM)
  a.xstructure2vasp();
  if (vasp5) {
    a.is_vasp4_poscar_format = false;
    a.is_vasp5_poscar_format = true;
  } // CO20210119
  //  cerr << a.title << endl;
  return a;
}

xstructure input2ITCxstr(istream& input) { // CO20220613
  xstructure a(input, IOAFLOW_AUTO);
  a.xstructure2itc();
  return a;
}

xstructure input2ELKxstr(istream& input) { // DX20200313
  xstructure a(input, IOAFLOW_AUTO);
  a.xstructure2elk();
  return a;
}

xstructure input2ATATxstr(istream& input) { // SD20220123
  xstructure a(input, IOAFLOW_AUTO);
  a.xstructure2atat();
  return a;
}

// **************************************************************************
// r_lattice
// ***************************************************************************
xvector<double> r_lattice(const xstructure& str, const int& l) {
  xvector<double> rrr(3);
  int i;
  int j;
  int k;
  l2ijk(str, l, i, j, k);
  rrr = ((double) i) * str.lattice(1) + ((double) j) * str.lattice(2) + ((double) k) * str.lattice(3);
  return rrr;
}

xvector<double> r_lattice(const xstructure& str, const int& i, const int& j, const int& k) {
  xvector<double> rrr(3);
  rrr = ((double) i) * str.lattice(1) + ((double) j) * str.lattice(2) + ((double) k) * str.lattice(3);
  return rrr;
}

xvector<double> r_lattice(const xstructure& str, const xvector<int>& ijk) {
  xvector<double> rrr(3);
  rrr = ((double) ijk(1)) * str.lattice(1) + ((double) ijk(2)) * str.lattice(2) + ((double) ijk(3)) * str.lattice(3);
  return rrr;
}

// **************************************************************************
// Function GetUnitCellRep
// **************************************************************************
// Dane Morgan / Stefano Curtarolo
// This function calculated the representation of a point in terms of a position in cell 0 and its true unit cell.
// Output position in cell 0 is given with same coordinate type (Cart or Direct) as input position.
// Note that there is always an ambiguity of being in a given cell with position 0, or one cell over with position 1.  This
// is broken by forcing all cell positions to be at 0 if they are within TOL of 1.  This gives consistent cell image locations.
// CO20170717 - Compare this with PBC() function, this function brings in to 0th cell, PBC brings in to cell between -0.5 and 0.5
// PBC() is good for minimizing overall fpos
void GetUnitCellRep(const xvector<double>& ppos, xvector<double>& p_cell0, xvector<int>& ijk, const xmatrix<double>& lattice, const bool coord_flag) {
  const double TOL = 1e-11;
  const xvector<double> cpos(3);
  xvector<double> fpos(3);
  // p_cell0=xvector<double> (3);
  // ijk=xvector<int> (3);

  clear(ijk);
  clear(p_cell0);

  if (coord_flag == true) { // ppos is in cartesian coords
    fpos = C2F(lattice, ppos);
  }

  for (int ic = 1; ic <= 3; ic++) {
    int s;
    s = SignNoZero(fpos(ic));
    if (s > 0) {
      ijk(ic) = int(fpos(ic));
    } else {
      ijk(ic) = int(fpos(ic)) - 1;
    }
    p_cell0(ic) = fpos(ic) - ijk(ic);
    // Boundary Correction.  p_cell0 is >=0 and <1.
    // If p_cell0 is within TOL of 1
    // then set it to 0 and shift ijk up by 1.
    if (std::abs(1 - p_cell0(ic)) < TOL) {
      p_cell0(ic) = 0;
      ijk(ic) = ijk(ic) + 1;
    }
  }
  if (coord_flag == true) { // ppos is in cartesian coords
    p_cell0 = F2C(lattice, p_cell0);
  }
}

// **************************************************************************
// Function COMPARE_GetNeighData
// **************************************************************************
// For sort algorithm in GetNeighData
class compare_GetNeighData {
public:
  int operator()(const _atom& a, const _atom& b) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << "operator()(const _atom& a, const _atom& b)" << endl;
    }
    const double tol = 1e-15;
    // Sort by distance
    if (aurostd::isequal(GetDistFromOrigin(a), GetDistFromOrigin(b), tol)) {
      if (LDEBUG) {
        cerr << "operator()(const _atom& a, const _atom& b)  isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)" << endl;
      }
      // Sort by unit cell values
      if (a.name == b.name) {
        if (LDEBUG) {
          cerr << "operator()(const _atom& a, const _atom& b)  a.name==b.name" << endl;
        }
        const xvector<int> ijka = a.ijk;
        const xvector<int> ijkb = b.ijk;
        const int va = 100 * ijka(1) + 10 * ijka(2) + 1 * ijka(3);
        const int vb = 100 * ijkb(1) + 10 * ijkb(2) + 1 * ijkb(3);
        return va < vb;
      } else {
        if (LDEBUG) {
          cerr << "operator()(const _atom& a, const _atom& b)  a.name!=b.name" << endl;
        }
        return a.name < b.name;
      }
    } else {
      if (LDEBUG) {
        cerr << "operator()(const _atom& a, const _atom& b)  !isequal(GetDistFromOrigin(a),GetDistFromOrigin(b),tol)" << endl;
      }
      return GetDistFromOrigin(a) < GetDistFromOrigin(b);
    }
    // Sort by name
    // if(a.name==b.name) {
    //  return GetDistFromOrigin(a)<GetDistFromOrigin(b);
    // }
    // else {
    //  return a.name<b.name;
    // }
  }
};

// RF20200831 - checkStructure - START
//  **************************************************************************
//  checkStructure
//  **************************************************************************
//  rescale structure to 1 and check whether e.g. species and atoms are present
void xstructure::checkStructure() {
  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;
  (*this).ReScale(1.0); // rescales scaling factor in second line of POSCAR to 1, needed for correct distances
  // throw some general information such as input structure
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << endl << "INPUT STRUCTURE:" << endl;
    cerr << __AFLOW_FUNC__ << (*this) << endl;
  }
  // check whether the species vector is populated, otherwise throw error
  if ((*this).species.empty()) {
    message << " BAD NEWS: It seems there are no species in the structure. Please adjust the structure and rerun.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
  }
  // check whether there are any atoms in the structure
  if ((*this).atoms.empty()) {
    message << " BAD NEWS: It seems there are no atoms in the structure. Please adjust the structure and rerun.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
  }
  // if species of atoms are not known like in VASP4 format, throw error
  for (size_t k = 0, ksize = (*this).atoms.size(); k < ksize; k++) {
    if ((*this).atoms[k].cleanname.empty()) {
      message << " BAD NEWS: It seems you are providing a structure without complete species information as input. This implementation requires a structure with the species information included. For a VASP4 "
                 "POSCAR, the species must be written on the right side next to the coordinates for each atom. Please adjust the structure and rerun.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
  }
}

// RF20200831 - checkStructure - END

// **************************************************************************
// Function GetNeighbors
// **************************************************************************
// rewrite of GetNeighData()
// atoms_cell is the atoms for which neighbors are found (could be iatoms), hence it determines the sizes of i_neighbors and distances
// i_neighbors are indices to grid_atoms (xstructure attribute)
// CO20200912
void xstructure::GetNeighbors(deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, double rmin, bool prim, bool unique_only) {
  deque<_atom> atoms_cell;
  return GetNeighbors(atoms_cell, i_neighbors, distances, rmin, prim, unique_only);
}

void xstructure::GetNeighbors(deque<_atom>& atoms_cell, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, double rmin, bool prim, bool unique_only) {
  const double rmax = RadiusSphereLattice((*this).scale * (*this).lattice);
  return GetNeighbors(atoms_cell, i_neighbors, distances, rmax, rmin, prim, unique_only);
}

void xstructure::GetNeighbors(deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, double rmax, double rmin, bool prim, bool unique_only) {
  deque<_atom> atoms_cell;
  return GetNeighbors(atoms_cell, i_neighbors, distances, rmax, rmin, prim, unique_only);
}

void xstructure::GetNeighbors(deque<_atom>& atoms_cell, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, double rmax, double rmin, bool prim, bool unique_only) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  uint i = 0;
  uint k = 0;

  if (prim) {
    GetPrimitive();
  }
  ReScale(1.0);
  // get atoms_cell
  atoms_cell.clear(); // clear
  for (i = 0; i < i_neighbors.size(); i++) {
    i_neighbors[i].clear();
  }
  i_neighbors.clear(); // clear
  for (i = 0; i < distances.size(); i++) {
    distances[i].clear();
  }
  distances.clear(); // clear
  vector<uint> atomscell2atoms_mapping;
  if (unique_only) {
    if (iatoms_calculated == false) {
      CalculateSymmetry();
    }
    for (i = 0; i < iatoms.size(); i++) {
      atoms_cell.push_back(atoms[iatoms[i][0]]);
      atomscell2atoms_mapping.push_back(iatoms[i][0]);
      i_neighbors.emplace_back(0);
      distances.emplace_back(0);
    }
  } else {
    for (size_t i = 0; i < atoms.size(); i++) {
      atoms_cell.push_back(atoms[i]);
      atomscell2atoms_mapping.push_back(i);
      i_neighbors.emplace_back(0);
      distances.emplace_back(0);
    }
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " atoms_cell.size()=" << atoms_cell.size() << endl;
    cerr << __AFLOW_FUNC__ << " atomscell2atoms_mapping=" << aurostd::joinWDelimiter(atomscell2atoms_mapping, ",") << endl;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " GenerateGridAtoms(): starting" << endl;
  }
  GenerateGridAtoms(rmax);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " GenerateGridAtoms(): done" << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " grid_atoms_number=" << grid_atoms_number << endl;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " generating neighbors" << endl;
  }
  double dist = 0.0;
  uint atom1 = 0;
  uint atom2 = 0;
  for (i = 0; i < atoms_cell.size(); i++) {
    atom1 = grid_atoms_pc2scMap[atomscell2atoms_mapping[i]];
    for (atom2 = 0; atom2 < (uint) grid_atoms_number; atom2++) {
      if (atom1 == atom2) {
        continue;
      } // skip self
      dist = AtomDist(grid_atoms[atom1], grid_atoms[atom2]); // distance
      if (dist > rmin && dist <= rmax) {
        i_neighbors[i].push_back(atom2);
        distances[i].push_back(dist);
      }
    }
  }

  // now sort
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " sorting" << endl;
  }
  for (k = 0; k < i_neighbors.size(); k++) {
    aurostd::sort(distances[k], i_neighbors[k]);
  }

  if (LDEBUG) {
    for (k = 0; k < i_neighbors.size(); k++) {
      cerr << __AFLOW_FUNC__ << " ATOMS_CELL[k=" << k << "]: " << atoms_cell[k].name << endl;
      for (i = 0; i < i_neighbors[k].size(); i++) {
        cerr << "  neighbor[i=" << i << "]: " << grid_atoms[i_neighbors[k][i]].name << " dist=" << distances[k][i] << endl;
      }
    }
  }
}

void GetNeighbors(const xstructure& xstr_in, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, double rmin, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetNeighbors(i_neighbors, distances, rmin, prim, unique_only);
}

void GetNeighbors(const xstructure& xstr_in, deque<_atom>& atoms_cell, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, double rmin, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetNeighbors(atoms_cell, i_neighbors, distances, rmin, prim, unique_only);
}

void GetNeighbors(const xstructure& xstr_in, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, double rmax, double rmin, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetNeighbors(i_neighbors, distances, rmax, rmin, prim, unique_only);
}

void GetNeighbors(const xstructure& xstr_in, deque<_atom>& atoms_cell, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, double rmax, double rmin, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetNeighbors(atoms_cell, i_neighbors, distances, rmax, rmin, prim, unique_only);
}

// **************************************************************************
// Function GetCoordinations
// **************************************************************************
// CO20200912
void xstructure::GetCoordinations(deque<deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only) {
  deque<deque<uint>> i_neighbors;
  deque<deque<double>> distances;
  return GetCoordinations(i_neighbors, distances, coordinations, tol, rmin, prim, unique_only);
}

void xstructure::GetCoordinations(deque<_atom>& atoms_cell, deque<deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only) {
  deque<deque<uint>> i_neighbors;
  deque<deque<double>> distances;
  return GetCoordinations(atoms_cell, i_neighbors, distances, coordinations, tol, rmin, prim, unique_only);
}

void xstructure::GetCoordinations(deque<deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only) {
  deque<deque<uint>> i_neighbors;
  deque<deque<double>> distances;
  return GetCoordinations(i_neighbors, distances, coordinations, rmax, rmin, tol, prim, unique_only);
}

void xstructure::GetCoordinations(deque<_atom>& atoms_cell, deque<deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only) {
  deque<deque<uint>> i_neighbors;
  deque<deque<double>> distances;
  return GetCoordinations(atoms_cell, i_neighbors, distances, coordinations, rmax, rmin, tol, prim, unique_only);
}

void xstructure::GetCoordinations(deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, deque<deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only) {
  deque<_atom> atoms_cell;
  return GetCoordinations(atoms_cell, i_neighbors, distances, coordinations, tol, rmin, prim, unique_only);
}

void xstructure::GetCoordinations(deque<_atom>& atoms_cell, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, deque<deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only) {
  const double rmax = RadiusSphereLattice((*this).scale * (*this).lattice);
  return GetCoordinations(atoms_cell, i_neighbors, distances, coordinations, rmax, rmin, tol, prim, unique_only);
}

void xstructure::GetCoordinations(deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, deque<deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only) {
  deque<_atom> atoms_cell;
  return GetCoordinations(atoms_cell, i_neighbors, distances, coordinations, rmax, rmin, tol, prim, unique_only);
}

void xstructure::GetCoordinations(deque<_atom>& atoms_cell, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, deque<deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  GetNeighbors(atoms_cell, i_neighbors, distances, rmax, rmin, prim, unique_only);

  uint i = 0;
  uint j = 0;
  uint k = 0;
  for (i = 0; i < coordinations.size(); i++) {
    coordinations[i].clear();
  }
  coordinations.clear(); // clear

  for (k = 0; k < i_neighbors.size(); k++) {
    coordinations.emplace_back(0);
    j = 0; // index to start at
    coordinations.back().push_back(0);
    for (i = 0; i < i_neighbors[k].size(); i++) {
      if (aurostd::isequal(distances[k][j], distances[k][i], tol)) {
        coordinations.back().back() += 1;
      } else {
        j = i;
        coordinations.back().push_back(1);
      }
    }
  }

  if (LDEBUG) {
    for (k = 0; k < i_neighbors.size(); k++) {
      cerr << __AFLOW_FUNC__ << " ATOMS_CELL[k=" << k << "]: " << atoms_cell[k].name << endl;
      for (i = 0; i < coordinations[k].size(); i++) {
        cerr << "  coordination[i=" << i << "]: " << coordinations[k][i] << endl;
      }
    }
  }
}

void GetCoordinations(const xstructure& xstr_in, deque<deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(coordinations, tol, rmin, prim, unique_only);
}

void GetCoordinations(const xstructure& xstr_in, deque<_atom>& atoms_cell, deque<deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(atoms_cell, coordinations, tol, rmin, prim, unique_only);
}

void GetCoordinations(const xstructure& xstr_in, deque<deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(coordinations, rmax, rmin, tol, prim, unique_only);
}

void GetCoordinations(const xstructure& xstr_in, deque<_atom>& atoms_cell, deque<deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(atoms_cell, coordinations, rmax, rmin, tol, prim, unique_only);
}

void GetCoordinations(const xstructure& xstr_in, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, deque<deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(i_neighbors, distances, coordinations, tol, rmin, prim, unique_only);
}

void GetCoordinations(const xstructure& xstr_in, deque<_atom>& atoms_cell, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, deque<deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(atoms_cell, i_neighbors, distances, coordinations, rmin, tol, prim, unique_only);
}

void GetCoordinations(const xstructure& xstr_in, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, deque<deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(i_neighbors, distances, coordinations, rmax, rmin, tol, prim, unique_only);
}

void GetCoordinations(const xstructure& xstr_in, deque<_atom>& atoms_cell, deque<deque<uint>>& i_neighbors, deque<deque<double>>& distances, deque<deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only) {
  xstructure xstr(xstr_in);
  return xstr.GetCoordinations(atoms_cell, i_neighbors, distances, coordinations, rmax, rmin, tol, prim, unique_only);
}

// **************************************************************************
// Function GetNeighData
// **************************************************************************
// This function collects all the neighbor data between rmin and rmax and stores
// it for each atom in a vector of atom objects in order of increasing distance.

void xstructure::GetNeighData(const double rmax, deque<deque<_atom>>& neigh_mat, const double rmin) const { // CO20220623
  deque<_atom> atoms_cell;
  return GetNeighData(atoms_cell, rmax, neigh_mat, rmin);
}

void xstructure::GetNeighData(deque<_atom>& atoms_cell, const double rmax, deque<deque<_atom>>& neigh_mat, const double rmin) const {
  // CO20220623 - rewritten using GetNeighbors() (MUCH faster)
  const bool LDEBUG = (false || XHOST.DEBUG);
  xstructure xstr_tmp(*this); // make copy so we can modify xstr
  deque<deque<uint>> i_neighbors;
  deque<deque<double>> distances;
  const bool prim = true;
  const bool unique_only = false;
  xstr_tmp.GetNeighbors(atoms_cell, i_neighbors, distances, rmax, rmin, prim, unique_only);
  // now build neigh_mat
  uint i = 0;
  uint j = 0;
  for (i = 0; i < neigh_mat.size(); i++) {
    neigh_mat[i].clear();
  }
  neigh_mat.clear(); // clear
  for (i = 0; i < i_neighbors.size(); i++) {
    neigh_mat.emplace_back(0);
    neigh_mat[i].push_back(atoms_cell[i]); // first atom is self
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " neigh_mat[i=" << i << "][j=0].cpos=" << neigh_mat[i].back().cpos << endl;
    }
    for (j = 0; j < i_neighbors[i].size(); j++) {
      neigh_mat[i].push_back(xstr_tmp.grid_atoms[i_neighbors[i][j]]);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " neigh_mat[i=" << i << "][j=" << j << "].cpos=" << neigh_mat[i].back().cpos << "; dist=" << distances[i][j] << endl;
      }
    }
  }
}

void xstructure::GetNeighData_20220101(const deque<_atom>& in_atom_vec, const double rmin, const double rmax, deque<deque<_atom>>& neigh_mat) const {
  // CO20220623 - AVOID USING THIS FUNCTION, use the one above
  // this one relies on an input of in_atom_vec which is usually the atoms of the structure
  // there is only one exception: GetGoodShellPoints() in aflow_pflow_print.cpp; this function should be rewritten in the future
  // considering the name of this function, it should only operate on the atoms of the structure
  // no point of having two sources of truth, it's confusing
  // the rewrite for an input of atoms should avoid ad-hoc, a posteriori mappings introducing noise in positions: GetUnitCellRep() and ConvertAtomToLat()
  const double epsilon = 1.0e-3; // sometimes you have wrong images due to roundoff
  const deque<_atom> neigh_vec;
  // Get data from str.
  xstructure sstr(*this);
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr.ReScale(1.0);
  xmatrix<double> original_lattice(3, 3);
  original_lattice = sstr.lattice;
  // Use Niggli cell to avoid missing neighbors in angular cells.
  // This involves converting your atoms into the 0th Niggli cell,
  // keeping track of the unit cell you came from, and then shifting back
  // at the end.  This shifting must be done for both atom_vec (at the
  // beginning) and neigh_mat (at the end).
  // sstr.NiggliUnitCellForm();
  sstr.MinkowskiBasisReduction();
  //    cerr << sstr << endl;

  xmatrix<double> lat(3, 3);
  lat = sstr.lattice;

  //  Convert input atom_vec to Niggli lattice
  deque<_atom> atom_vec = in_atom_vec;
  const xmatrix<int> atom_shifts((int) atom_vec.size() - 1, 3, 0, 1);
  for (size_t ia = 0; ia < atom_vec.size(); ia++) {
    const _atom a = atom_vec[ia];
    xvector<double> p_cell0(3);
    xvector<int> atom_shifts_temp(3);
    GetUnitCellRep(a.cpos, p_cell0, atom_shifts_temp, lat, true);
    atom_vec[ia].cpos = p_cell0;
    atom_vec[ia] = ConvertAtomToLat(atom_vec[ia], lat);
    atom_shifts(ia, 1) = atom_shifts_temp(1);
    atom_shifts(ia, 2) = atom_shifts_temp(2);
    atom_shifts(ia, 3) = atom_shifts_temp(3);
  }

  // Create vector of all atoms in all periodic images needed so that for atoms in the unit cell every possible neighbor within rmax
  // is included.  This can be done by looping over all atoms in all unit cells such that if the unit cells are given by
  // n0*v0,n1*v1,n2*v2, then ni_max*vi>rmax for all i=0->2.  This is not rigorous but seems to work quite well.  It failed for very
  // angular cells, but by using the Niggli reduced cell this pitfall seems to be avoided (but no proof, so be careful).

  int imax;
  int jmax;
  int kmax;
  // Find imax
  // (algorithm 1) approximate
  imax = (int) max(rmax / aurostd::modulus(lat(1)), rmax / aurostd::modulus(lat(2)), rmax / aurostd::modulus(lat(3))) + 2;
  jmax = imax;
  kmax = imax;
  // (algorithm 2) exact, the requirement is that atoms must in incell.
  // we should find nlat1,nlat2,nlat3
  // where nlat1 is the lat1 component that is perpendicular to the plane made by lat2 & lat3.
  // then imax=ceil(rmax/nlat1)
  // This algorithm is implemented in function LatticeDimensionSphere
  xvector<int> dims(3);
  dims = LatticeDimensionSphere(lat, rmax);
  imax = dims(1);
  jmax = dims(2);
  kmax = dims(3);
  const xvector<int> ijk(3);
  deque<_atom> all_atom_vec;
  // latt maybe a rotated version of POSCAR, so need to incell-ized the fpos and cpos
  // sstr.BringInCell(); // with roundoff
  sstr.BringInCell(); // DX+RF20200121 - negative input is no longer accepted // no roundoff
  // DX+RF20200121 - negative input is no longer accepted : sstr.BringInCell(-1.0); // no roundoff
  for (ijk(1) = -imax; ijk(1) <= imax; ijk(1)++) {
    for (ijk(2) = -jmax; ijk(2) <= jmax; ijk(2)++) {
      for (ijk(3) = -kmax; ijk(3) <= kmax; ijk(3)++) {
        const xvector<double> ctpos(3);
        xvector<double> ftpos(3);
        for (size_t iat = 0; iat < sstr.atoms.size(); iat++) {
          _atom a;
          a = sstr.atoms[iat];
          a.name = sstr.atoms[iat].name;
          a.basis = iat; //[CO20200130 - number->basis]a.number=iat;
          a.ijk = ijk;
          for (int ic = 1; ic <= 3; ic++) {
            ctpos(ic) = sstr.atoms[iat].cpos(ic) + ijk(1) * lat(1, ic) + ijk(2) * lat(2, ic) + ijk(3) * lat(3, ic);
          }
          ftpos = C2F(lat, ctpos);
          a.cpos = ctpos;
          a.fpos = ftpos;
          a.type = sstr.atoms[iat].type;
          all_atom_vec.push_back(a);
        } // iat
      } // ijk
    } // ijk
  } // ijk

  // Now build neighbors list for each atom on atom list.  Each
  // neighbor list will be row in the neigh_mat matrix.
  for (size_t ia1 = 0; ia1 < atom_vec.size(); ia1++) {
    xvector<double> corigin(3);
    _atom at = atom_vec[ia1];
    corigin = at.cpos;
    at.corigin = corigin;
    deque<_atom> neigh_vec;
    neigh_vec.push_back(at);
    for (size_t ia2 = 0; ia2 < all_atom_vec.size(); ia2++) {
      const double dist = AtomDist(at, all_atom_vec[ia2]);
      if (dist <= rmax && dist >= rmin && dist >= epsilon) {
        all_atom_vec[ia2].corigin = corigin;
        neigh_vec.push_back(all_atom_vec[ia2]);
      } // if
    } // ia2
    sort(neigh_vec.begin() + 1, neigh_vec.end(), compare_GetNeighData());
    neigh_mat.push_back(neigh_vec);
  } // ia1

  //  Convert neigh_mat from Niggli cell to original lattice
  for (size_t ia = 0; ia < neigh_mat.size(); ia++) {
    for (size_t ja = 0; ja < neigh_mat[ia].size(); ja++) {
      xvector<double> fpos(3);
      xvector<double> cpos(3);
      _atom a;
      a = neigh_mat[ia][ja];
      fpos = a.fpos;
      for (int ic = 1; ic <= 3; ic++) {
        fpos(ic) = fpos(ic) + atom_shifts(ia, ic);
      }
      cpos = F2C(lat, fpos);
      a.cpos = cpos;
      neigh_mat[ia][ja] = ConvertAtomToLat(a, original_lattice);
    }
  }
}

// **************************************************************************
// Function GetBasisTransformation //DX20201015
// **************************************************************************
xmatrix<double> GetBasisTransformation(const xmatrix<double>& lattice_original, const xmatrix<double>& lattice_new) {
  return lattice_new * inverse(lattice_original);
}

// **************************************************************************
// Function GetBasisTransformationInternalTranslations //DX20201124
// **************************************************************************
vector<xvector<double>> GetBasisTransformationInternalTranslations(const xmatrix<double>& basis_transformation) {
  // Given a basis transformation matrix, determine any internal lattice
  // translation(s). This is necessary if the basis transformation increases
  // the volume of the cell, otherwise, there are no internal translations
  // (return immediately).
  // Another way to think of this: if you expand your lattice/cell, this
  // function finds all the lattice points in the new cell

  const bool LDEBUG = (false || XHOST.DEBUG);

  vector<xvector<double>> translations;

  // ---------------------------------------------------------------------------
  // check if the basis transformation makes the cell larger and find
  // corresponding internal translations
  // DX20210520 - DO NOT EXCLUDE BASED ON DETERMINANT
  // It is possible to stretch in one direction and compress in another and still
  // get determinant=1 (e.g., POCC structures)
  // DX20210520 [OBSOELTE] if(cell_volume_change-1.0>_AUROSTD_XSCALAR_TOLERANCE_INTEGER_){}

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " cell size increases. Finding internal translations." << endl;
  }

  // ---------------------------------------------------------------------------
  // get inverse matrix (Q)
  const xmatrix<double> inverse_transform = aurostd::inverse(basis_transformation);

  // ---------------------------------------------------------------------------
  // to get translations take the "larger cell" in fractional coordinates
  // and perform the inverse operation (Q) to see how small it gets,
  // then these are the internal translations
  const xmatrix<double> lattice_frac = aurostd::eye<double>(3, 3);
  const xmatrix<double> lattice_shrink = inverse_transform * lattice_frac;

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " shrunken lattice: " << lattice_shrink << endl;
  }

  // ---------------------------------------------------------------------------
  // Now that we have the shortest internal translations from lattice shrink
  // (forms a basis), we need to find all the internal translations inside this
  // cell via linear combinations of this basis.
  // To determine how many combinations we need (i.e. how far to expand), we can
  // use LatticeDimensionSphere(). Since lattice_shrink is in fractional
  // coordinates, we need to find the necessary dimensions in each direction
  // to fill the cell (i.e., the unit box). //DX20210111
  const xvector<int> dims = LatticeDimensionSphere(lattice_shrink, 1.0);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " number of times to apply each internal translation: " << dims[1] << "," << dims[2] << "," << dims[3] << endl;
  }

  // ---------------------------------------------------------------------------
  // create all linear combinations of translations, filter out duplicates later
  const xvector<double> a_vec = lattice_shrink(1);
  const xvector<double> b_vec = lattice_shrink(2);
  const xvector<double> c_vec = lattice_shrink(3);
  xvector<double> a_vec_scaled;
  xvector<double> b_vec_scaled;
  xvector<double> c_vec_scaled;
  for (int a = 0; a <= dims[1]; a++) { // DX20210506 - need <=
    a_vec_scaled = (double) a * a_vec;
    translations.push_back(a_vec_scaled);
    for (int b = 0; b <= dims[2]; b++) { // DX20210506 - need <=
      b_vec_scaled = (double) b * b_vec;
      translations.push_back(b_vec_scaled);
      translations.push_back(a_vec_scaled + b_vec_scaled);
      for (int c = 0; c <= dims[3]; c++) { // DX20210506 - need <=
        c_vec_scaled = (double) c * c_vec;
        translations.push_back(c_vec_scaled);
        translations.push_back(a_vec_scaled + c_vec_scaled);
        translations.push_back(b_vec_scaled + c_vec_scaled);
        translations.push_back(a_vec_scaled + b_vec_scaled + c_vec_scaled);
      }
    }
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " # translations:" << translations.size() << endl;
    for (size_t t = 0; t < translations.size(); t++) {
      cerr << __AFLOW_FUNC__ << " translations:" << translations[t] << endl;
    }
  }

  // ---------------------------------------------------------------------------
  // filter out unique translations
  vector<xvector<double>> unique_translations;
  bool unique = true;
  for (size_t t = 0; t < translations.size(); t++) {
    const xvector<double> translation_incell = BringInCell(translations[t]);
    unique = true;
    for (size_t u = 0; u < unique_translations.size() && unique; u++) {
      unique = !(aurostd::isequal(translation_incell, unique_translations[u]));
    }
    if (unique) {
      unique_translations.push_back(translation_incell);
    }
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " # unique_translations:" << unique_translations.size() << endl;
    for (size_t t = 0; t < unique_translations.size(); t++) {
      cerr << __AFLOW_FUNC__ << " unique_translations:" << unique_translations[t] << endl;
    }
  }
  translations = unique_translations;

  // ---------------------------------------------------------------------------
  // if the cell size remains the same or shrinks, no internal translations
  // DX20210520 [OBSOELTE] else{}
  // DX20210520 [OBSOELTE]  // use null vector
  // DX20210520 [OBSOELTE]  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " cell size remains the same or reduced. No internal translations." << endl; }
  // DX20210520 [OBSOELTE]  xvector<double> zero_xvector;
  // DX20210520 [OBSOELTE]  translations.push_back(zero_xvector);
  return translations;
}

// **************************************************************************
// Function GetRotation //DX20201015
// **************************************************************************
xmatrix<double> GetRotation(const xmatrix<double>& lattice_original, const xmatrix<double>& lattice_new) {
  return aurostd::inverse(lattice_original) * lattice_new;
}

// **************************************************************************
// Function ChangeBasis() //DX20201015
// **************************************************************************
// Convert a structure (lattice and atom positions) into a new representation
// based on the input transformation matrix.
// The transformation matrix is generally NOT a unitary transformation -
// it can change the volume of the cell - otherwise it would be a rotation
// (use Rotate() instead).
// The procedure is generalized for transformations that enlarge (supercell)
// or reduce (primitivize) the structure.
// Enlarging the cell: search for unique internal translations based on
// transformation matrix.
// Reducing the cell: remove duplicate atom positions.
// Example transformation matrix (4x1x1 supercell expansion):
//   -4.0000e+00  0.0000e+00  0.0000e+00
//   -1.0000e+00  0.0000e+00  1.0000e+00
//   -1.0000e+00  1.0000e+00  0.0000e+00

// ---------------------------------------------------------------------------
// returns new xstructure (makes a copy)
xstructure ChangeBasis(const xstructure& xstr, const xmatrix<double>& transformation_matrix) {
  xstructure xstr_transformed = xstr;
  xstr_transformed.ChangeBasis(transformation_matrix);
  return xstr_transformed;
}

// ---------------------------------------------------------------------------
// modifies in-place (efficient)
void xstructure::ChangeBasis(const xmatrix<double>& transformation_matrix) {
  // if the transformation matrix is the identity, don't do anything
  if (aurostd::isidentity(transformation_matrix)) {
    return;
  }

  const bool LDEBUG = (false || XHOST.DEBUG);
  stringstream message;

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " structure BEFORE basis transformation:" << endl;
    cerr << (*this) << endl;
  }

  const size_t natoms_orig = (*this).atoms.size();
  uint natoms_transformed = 0;
  bool is_integer_multiple_transformation = true;
  if ((*this).dist_nn_min == AUROSTD_NAN || (*this).dist_nn_min == AUROSTD_MAX_DOUBLE) {
    (*this).dist_nn_min = SYM::minimumDistance((*this));
  }

  // if the cell is being reduced in a symmetry routine, it is best to set the same-atom tolerance to
  // the minimum interatomic distance minus the sym_eps (e.g., dist_nn_min-sym_eps), since this is the "new resolution"
  // of the atom postions; 90% of dist_nn_min is used so that the nearest-neighbor is not removed erroneously //DX20210716
  double tol = 0.9 * (*this).dist_nn_min;
  // DX20210316 - 0.1 is not robust, used fraction of min_dist //DX20210716 - default tol of 90% of dist_nn_min should be sufficient
  if ((*this).sym_eps != AUROSTD_NAN && (*this).sym_eps != AUROSTD_MAX_DOUBLE) {
    tol -= (*this).sym_eps;
  }
  // DX20210716 - if sym_eps is valid, then subtract it from the tol, since this is the new "resolution" for the system

  // ---------------------------------------------------------------------------
  // transform the lattice
  const xmatrix<double> lattice_orig = (*this).lattice;
  (*this).lattice = transformation_matrix * (*this).lattice;
  (*this).FixLattices();

  // ---------------------------------------------------------------------------
  // get internal translations from basis transformation (i.e. transforming
  // to larger cells)
  vector<xvector<double>> translations = GetBasisTransformationInternalTranslations(transformation_matrix);

  // ---------------------------------------------------------------------------
  // transform the atom positions
  deque<_atom> atom_basis;
  _atom atom_tmp;
  const xmatrix<double> forig2fnew = inverse(trasp(transformation_matrix));
  // Q*pos , but need to transpose Q for AFLOW xmatrix convention
  const xmatrix<double> f2c = trasp((*this).lattice); // Q*pos , but need to transpose Q for AFLOW xmatrix convention
  for (size_t i = 0; i < (*this).atoms.size(); i++) {
    for (size_t t = 0; t < translations.size(); t++) {
      atom_tmp = (*this).atoms[i];
      atom_tmp.fpos = forig2fnew * ((*this).atoms[i].fpos);
      atom_tmp.fpos = ::BringInCell(atom_tmp.fpos + translations[t]);
      atom_tmp.cpos = f2c * atom_tmp.fpos;
      if (!SYM::MapAtom(atom_basis, atom_tmp, true, (*this).lattice, false, tol)) {
        // DX20210324 - to account for duplicate translations
        atom_basis.push_back(atom_tmp);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // calculate change in basis transformation determinant
  // (shift to 1 to easily see if reduces or expands)
  const double basis_transformation_det_change = aurostd::abs(aurostd::det(transformation_matrix)) - 1.0;

  // ---------------------------------------------------------------------------
  // reduce the cell: remove any duplicate atoms
  // use _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ to be consistent with AUROSTD's
  // isinteger tolerance
  if (basis_transformation_det_change < -_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " removing duplicate atoms (cell has been reduced)." << endl;
    }

    const bool skew = false;
    const deque<_atom> new_basis = ::foldAtomsInCell(atom_basis, lattice_orig, (*this).lattice, skew, tol, false);
    // false: don't check atom mappings (slow) //DX20210118 - add global namespace
    atom_basis = new_basis;

    // check atom count
    natoms_transformed = atom_basis.size();
    is_integer_multiple_transformation = ((*this).num_each_type.size() == 1 || natoms_orig % natoms_transformed == 0);
    // DX20210316 - integer multiple does not apply to unaries
  }
  // ---------------------------------------------------------------------------
  // enlarge the cell: update the atom count information
  else if (basis_transformation_det_change > _AUROSTD_XSCALAR_TOLERANCE_INTEGER_) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " cell size has increased." << endl;
    }
    // check atom count
    natoms_transformed = atom_basis.size();
    is_integer_multiple_transformation = ((*this).num_each_type.size() == 1 || natoms_transformed % natoms_orig == 0);
    // DX20210316 - integer multiple does not apply to unaries
  }

  // ---------------------------------------------------------------------------
  // is integer multiple transformation (reduce or enlarge)
  if (!is_integer_multiple_transformation) {
    message << "Number of atoms is no longer an integer multiple with respect to the input structure"
            << " original: " << natoms_orig << " transformed: " << natoms_transformed << "; check the transformation matrix or same-atom tolerance.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  // ---------------------------------------------------------------------------
  // if the number of atoms changed (i.e., change in determinant is zero),
  // update the atom counts/order/types/etc.
  if (!aurostd::isequal(aurostd::abs(basis_transformation_det_change), _ZERO_TOL_, _AUROSTD_XSCALAR_TOLERANCE_INTEGER_)) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " updating atom count information." << endl;
    }
    std::stable_sort(atom_basis.begin(), atom_basis.end(), sortAtomsNames); // DX20210129
    (*this).ReplaceAtoms(atom_basis, false); // false: check_atom_overlap
  }
  // ---------------------------------------------------------------------------
  // if the transformation preserves the volume, one-to-one mappings
  // no need to update species/types/etc. (i.e., ReplaceAtoms() is not needed)
  else {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " cell size remains the same (updating atom positions)." << endl;
    }
    (*this).atoms = atom_basis;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " structure AFTER basis transformation:" << endl;
    cerr << (*this) << endl;
  }
}

// **************************************************************************
// Function TransformStructure() //DX20201125
// **************************************************************************
// ---------------------------------------------------------------------------
// returns new xstructure (makes a copy)
xstructure TransformStructure(const xstructure& xstr, const xmatrix<double>& transformation_matrix, const xmatrix<double>& rotation) {
  const xvector<double> origin_shift;
  return TransformStructure(xstr, transformation_matrix, rotation, origin_shift);
}

xstructure TransformStructure(const xstructure& xstr, const xmatrix<double>& transformation_matrix, const xmatrix<double>& rotation, const xvector<double>& origin_shift, bool is_shift_frac) {
  xstructure xstr_transformed = xstr;
  xstr_transformed.TransformStructure(transformation_matrix, rotation, origin_shift, is_shift_frac);
  return xstr_transformed;
}

// ---------------------------------------------------------------------------
// modifies in-place (efficient)
void xstructure::TransformStructure(const xmatrix<double>& transformation_matrix, const xmatrix<double>& rotation) {
  const xvector<double> origin_shift;
  (*this).TransformStructure(transformation_matrix, rotation, origin_shift);
}

void xstructure::TransformStructure(const xmatrix<double>& transformation_matrix, const xmatrix<double>& rotation, const xvector<double>& origin_shift, bool is_shift_frac) {
  const bool LDEBUG = (false || XHOST.DEBUG);

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " basis transformation: " << transformation_matrix << endl;
    cerr << __AFLOW_FUNC__ << " rotation (R): " << rotation << endl;
  }

  // ---------------------------------------------------------------------------
  // changed basis
  (*this).ChangeBasis(transformation_matrix);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " structure after CHANGING BASIS: " << (*this) << endl;
  }

  // ---------------------------------------------------------------------------
  // rotate
  (*this).Rotate(rotation);
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " structure after ROTATING: " << (*this) << endl;
  }

  // ---------------------------------------------------------------------------
  // rotate
  const bool coordinate_flag = (*this).coord_flag; // store original coordinate-type
  (*this).ShiftPos(origin_shift, is_shift_frac);
  (*this).coord_flag = coordinate_flag; // set back to original coordinate-type
  (*this).BringInCell(); // DX20210116
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " structure after shifting origin: " << (*this) << endl;
  }
}

// **************************************************************************
// GenerateXXXXXXXXXXx
// **************************************************************************
void xstructure::qm_clear() {
  qm_calculated = false;
  qm_scale = 0.0;
  qm_lattice.clear();
  qm_origin.clear();
  qm_klattice.clear();
  qm_f2c.clear();
  qm_c2f.clear();
  qm_origin.clear();
  qm_forces_write = false;
  qm_positions_write = false;
  qm_atoms.clear();
  qm_forces.clear();
  qm_positions.clear();
  for (size_t iat = 0; iat < atoms.size(); iat++) {
    qm_atoms.push_back(atoms[iat]); // just plug something, then I`ll clean
    qm_atoms.at(iat).cpos.clear();
    qm_atoms.at(iat).fpos.clear();
    qm_forces.push_back(atoms[iat].cpos); // just plug something, then I`ll clean
    qm_forces.at(iat).clear();
    qm_positions.push_back(atoms[iat].cpos); // just plug something, then I`ll clean
    qm_positions.at(iat).clear();
  }
  qm_E_cell = 0.0;
  qm_dE_cell = 0.0;
  qm_H_cell = 0.0;
  qm_PV_cell = 0.0;
  qm_mag_cell = 0.0;
  qm_P = 0.0;
  qm_E_atom = 0.0;
  qm_dE_atom = 0.0;
  qm_H_atom = 0.0;
  qm_PV_atom = 0.0;
  qm_mag_atom = 0.0;

  if (atoms.size() != qm_atoms.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_clear():", "[1] atoms.size()!=qm_atoms.size().", _VALUE_ERROR_);
  }
  if (atoms.size() != qm_forces.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_clear():", "[2] atoms.size()!=qm_forces.size().", _VALUE_ERROR_);
  }
  if (atoms.size() != qm_positions.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_clear():", "[3] atoms.size()!=qm_positions.size().", _VALUE_ERROR_);
  }
}

void xstructure::qm_recycle() {
  if (atoms.size() != qm_atoms.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_recycle():", "[1] atoms.size()!=qm_atoms.size().", _VALUE_ERROR_);
  }
  if (atoms.size() != qm_forces.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_recycle():", "[2] atoms.size()!=qm_forces.size().", _VALUE_ERROR_);
  }
  if (atoms.size() != qm_positions.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_recycle():", "[3] atoms.size()!=qm_positions.size().", _VALUE_ERROR_);
  }
  scale = qm_scale;
  lattice = qm_lattice;
  klattice = qm_klattice;
  f2c = qm_f2c;
  c2f = qm_c2f;
  origin = qm_origin;
  FixLattices();
  for (size_t i = 0; i < atoms.size(); i++) { // copy cpos/fpos
    atoms[i].cpos = qm_atoms.at(i).cpos; // get cpos
    atoms[i].fpos = qm_atoms.at(i).fpos; // get fpos
  }
  qm_clear();
}

void xstructure::qm_load(const string& Directory, const string& suffix, int iomode) {
  const double data_natoms = double(atoms.size());
  if (iomode != IOVASP_POSCAR) {
    throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", "Only IOVASP_POSCAR is supported.", _FILE_WRONG_FORMAT_);
  };
  if (iomode == IOVASP_POSCAR) {
    xOUTCAR outcar;
    if (aurostd::FileEmpty(Directory + "/OUTCAR" + suffix)) {
      throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", "Empty OUTCAR.", _FILE_CORRUPT_);
    } // PN+JJPR FIXED BUG

    outcar.GetPropertiesFile(Directory + "/OUTCAR" + suffix); // SD20221019 - don't perform natoms check yet

    if (std::abs(data_natoms - outcar.natoms) > 0.1) {
      stringstream message;
      message << "ERROR void xstructure::qm_load: data_natoms(" << data_natoms << ")!= (int) outcar.natoms(" << outcar.natoms << ") ..." << endl;
      message << "      Directory=" << Directory << endl;
      message << "      suffix=" << suffix << endl;
      message << "      iomode=" << iomode << endl;
      throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", message, _FILE_WRONG_FORMAT_);
      ;
    }

    //    cerr << atoms.size() << endl;
    qm_clear();
    if (atoms.size() != qm_atoms.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", "[1] atoms.size()!=qm_atoms.size().", _VALUE_ERROR_);
    }
    if (atoms.size() != qm_forces.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", "[2] atoms.size()!=qm_forces.size().", _VALUE_ERROR_);
    }
    if (atoms.size() != qm_positions.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", "[3] atoms.size()!=qm_positions.size().", _VALUE_ERROR_);
    }

    // NEW WITH xOUTCAR
    qm_forces.clear();
    for (size_t i = 0; i < outcar.vforces.size(); i++) {
      qm_forces.push_back(outcar.vforces[i]);
    }
    qm_positions.clear();
    for (size_t i = 0; i < outcar.vpositions_cartesian.size(); i++) {
      qm_positions.push_back(outcar.vpositions_cartesian[i]);
    }
    if (atoms.size() != qm_forces.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", "[4] atoms.size()!=qm_forces.size().", _VALUE_ERROR_);
    }
    if (atoms.size() != qm_positions.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", "[5] atoms.size()!=qm_positions.size().", _VALUE_ERROR_);
    }

    // NEW WITH xVASPRUNXML
    xVASPRUNXML vasprunxml;
    if (aurostd::FileEmpty(Directory + "/vasprun.xml" + suffix)) {
      throw aurostd::xerror(__AFLOW_FILE__, XPID + "xstructure::qm_load():", "Empty vasprun.xml.", _FILE_CORRUPT_);
    } // PN+JJPR FIXED BUG
    vasprunxml.GetPropertiesFile(Directory + "/vasprun.xml" + suffix);
    qm_forces.clear();
    for (size_t i = 0; i < vasprunxml.vforces.size(); i++) {
      qm_forces.push_back(vasprunxml.vforces[i]);
    }

    // OLD
    stringstream CONTCAR_stringstream;
    CONTCAR_stringstream.str(std::string());
    CONTCAR_stringstream << aurostd::file2string(Directory + "/CONTCAR" + suffix);

    qm_E_cell = outcar.energy_cell;
    qm_E_atom = outcar.energy_atom;
    qm_H_cell = outcar.enthalpy_cell;
    qm_H_atom = outcar.enthalpy_atom;
    qm_PV_cell = outcar.PV_cell;
    qm_PV_atom = outcar.PV_atom;
    qm_P = outcar.pressure;
    qm_mag_cell = outcar.mag_cell;
    qm_mag_atom = outcar.mag_atom;

    // CONTCAR OPERATIONS ---------------------------------------------------------------
    xstructure b(CONTCAR_stringstream); // LOAD it in
    qm_scale = b.scale;
    qm_lattice = b.lattice;
    qm_klattice = b.klattice;
    qm_f2c = b.f2c;
    qm_c2f = b.c2f;
    qm_origin = b.origin;
    for (size_t i = 0; i < atoms.size(); i++) { // copy cpos/fpos
      qm_atoms.at(i) = atoms[i];
      qm_atoms.at(i).cpos = b.atoms.at(i).cpos; // get cpos
      qm_atoms.at(i).fpos = b.atoms.at(i).fpos; // get fpos
    }
    qm_calculated = true;
  }
}

/// @brief create a x3d representation of a xstructure that can be converted in different file formats
/// @note not optimized for supercells yet
/// @author
/// @mod{HE,20250505,created}
aurostd::x3DWriter xstructure::render() const {
  aurostd::x3DWriter w;
  constexpr double atom_size_factor = 10 * 0.5; // 10: Ang to nm; 0.5: don't render atoms space filling
  w.prepareSceneLattice(lattice);
  w.tachyon_camera_orthographic = true;
  w.addLatticeBox(lattice, 0.075, "bm_grey", true);
  std::map<std::string, std::string> species_materials = w.addNamedColorSpreadMaterial(species);
  for (const auto& atom : atoms) {
    const double radius = GetAtomRadius(atom.atomic_number);

    // check if the atom is on a cell border
    std::bitset<3> direction_array; // defaults to all false
    for (int idx = 1; idx < 4; idx++) {
      if (std::fabs(atom.fpos[idx] - 1) < 1E-4 or std::fabs(atom.fpos[idx]) < 1E-4) {
        direction_array[idx - 1] = true;
      }
    }

    // add atom and all its images along the cell border
    xvector<double> new_fpos;
    switch (direction_array.count()) {
      case 0: w.addSphere(atom.cpos, atom_size_factor * radius, species_materials[atom.cleanname]); break;
      case 1: // cell wall
        for (const double da : {0.0, 1.0}) {
          if (direction_array[0]) {
            new_fpos = {da, atom.fpos[2], atom.fpos[3]};
          } else if (direction_array[1]) {
            new_fpos = {atom.fpos[1], da, atom.fpos[3]};
          } else if (direction_array[2]) {
            new_fpos = {atom.fpos[1], atom.fpos[2], da};
          }
        }
        break;
      case 2: // cell edge
        for (const double da : {0.0, 1.0}) {
          for (const double db : {0.0, 1.0}) {
            if (!direction_array[0]) {
              new_fpos = {atom.fpos[1], da, db};
            } else if (!direction_array[1]) {
              new_fpos = {da, atom.fpos[2], db};
            } else if (!direction_array[2]) {
              new_fpos = {da, db, atom.fpos[3]};
            }
            w.addSphere(f2c * new_fpos, atom_size_factor * radius, species_materials[atom.cleanname]);
          }
        }
        break;
      case 3: // cell corner
        for (const double dx : {0.0, 1.0}) {
          for (const double dy : {0.0, 1.0}) {
            for (const double dz : {0.0, 1.0}) {
              new_fpos = {dx, dy, dz};
              w.addSphere(f2c * new_fpos, atom_size_factor * radius, species_materials[atom.cleanname]);
            }
          }
        }
        break;
      default: break; // impossible to reached
    }
  }
  return w;
}

/// @brief convert a xstructure to json
/// @param xstr xstructure to convert
/// @return json object representing `xstr`
///
/// @author
/// @mod{DX,20170831,created}
/// @mod{HE,20240221,change to aurostd::JSON::object}
aurostd::JSON::object xstructure2json(const xstructure& xstr) {
  aurostd::JSON::object json(aurostd::JSON::object_types::DICTIONARY);

  if (!xstr.title.empty()) {
    json["title"] = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title);
  } else {
    if (PRINT_NULL_JSON) {
      json["title"] = nullptr;
    }
  }

  if (xstr.scale) {
    json["scale"] = xstr.scale;
  } else {
    if (PRINT_NULL_JSON) {
      json["scale"] = nullptr;
    }
  }

  if (xstr.lattice.rows) {
    json["lattice"] = xstr.lattice;
  } else {
    if (PRINT_NULL_JSON) {
      json["lattice"] = nullptr;
    }
  }

  if (!xstr.species.empty()) {
    vector<string> cleaned_species; // DX20190612 - cleaned species names
    for (size_t i = 0; i < xstr.species.size(); i++) {
      cleaned_species.push_back(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i]));
    } // DX20190612 - cleaned species names
    json["species"] = cleaned_species;
  } else {
    if (PRINT_NULL_JSON) {
      json["species"] = nullptr;
    }
  }

  if (!xstr.num_each_type.empty()) {
    json["number_each_type"] = xstr.num_each_type;
  } else {
    if (PRINT_NULL_JSON) {
      json["number_each_type"] = nullptr;
    }
  }

  if (xstr.coord_flag == _COORDS_FRACTIONAL_) {
    json["coordinates_type"] = "direct";
  } else if (xstr.coord_flag == _COORDS_CARTESIAN_) {
    json["coordinates_type"] = "Cartesian";
  }

  // ATOMS
  if (!xstr.atoms.empty()) {
    json["atoms"] = aurostd::JSON::object(aurostd::JSON::object_types::LIST);
    for (size_t i = 0; i < xstr.atoms.size(); i++) {
      json["atoms"].push_back(atom2json(xstr.atoms[i], xstr.coord_flag, xstr.partial_occupation_flag));
    }
  } else {
    if (PRINT_NULL_JSON) {
      json["atoms"] = nullptr;
    }
  }

  return json;
}

/// @brief convert an atom object to JSON
/// @param atom to convert
/// @return json object representing `atom`
///
/// @author
/// @mod{DX,20170831,created}
/// @mod{HE,20240221,change to aurostd::JSON::object}
aurostd::JSON::object atom2json(const _atom& atom, int coord_flag, int poccupation) {
  const bool roff = true; // round off
  const vector<string> vcontent_json;

  aurostd::JSON::object json(aurostd::JSON::object_types::DICTIONARY);

  if (!atom.name.empty()) {
    json["name"] = KBIN::VASP_PseudoPotential_CleanName(atom.name);
    // DX20190612 - added function to clean names
  } else {
    if (PRINT_NULL_JSON) {
      json["name"] = nullptr;
    }
  }

  if (coord_flag == _COORDS_FRACTIONAL_) {
    json["position"] = atom.fpos;
  } else if (coord_flag == _COORDS_CARTESIAN_) {
    json["position"] = atom.cpos;
  } else {
    if (PRINT_NULL_JSON) {
      json["position"] = nullptr;
    }
  }

  if (poccupation == true) {
    json["occupancy"] = atom.partial_occupation_value;
  } else if (poccupation == false) {
    json["occupancy"] = 1.0;
  } else {
    if (PRINT_NULL_JSON) {
      json["occupancy"] = nullptr;
    }
  }

  return json;
}

namespace soliquidy {
  /// @namespace soliquidy
  /// @brief tools to calculate Soliquidy
  ///
  /// @authors
  /// @mod{HE,20240530,created}
  ///
  /// @see
  /// @xlink{"Soliquidy: a geometric descriptor for phase transitions",https://materials.duke.edu/publications}
  /// @doi{10.1016/j.cag.2018.01.009,"Notions of optimal transport theory and how to implement them on a computer"}

  /// @brief execute soliquidy code based on comandline input
  /// @param vpflow option collection (SOLIQUIDY::AUID, SOLIQUIDY::WORKLIST, SOLIQUIDY::OUTPUT, SOLIQUIDY::X3D)
  /// @param input standard input
  void run_cmd(aurostd::xoption& vpflow, istream& input) {
    if (vpflow.flag("SOLIQUIDY::AUID")) {
      std::vector<std::string> auid_list;
      aurostd::string2tokens(vpflow.getattachedscheme("SOLIQUIDY::AUID"), auid_list, ",");
      cout << vpflow.getattachedscheme("SOLIQUIDY::OUTPUT") << endl;
      cout << CalculateAUID(auid_list, vpflow.getattachedscheme("SOLIQUIDY::OUTPUT"), vpflow.flag("SOLIQUIDY::X3D")).toString() << endl;
    } else if (vpflow.flag("SOLIQUIDY::WORKLIST")) {
      const std::string file_content = aurostd::file2string(vpflow.getattachedscheme("SOLIQUIDY::WORKLIST"));
      std::vector<std::string> worklist;
      aurostd::string2tokens(file_content, worklist, "\n");
      cout << CalculateWorklist(worklist, vpflow.getattachedscheme("SOLIQUIDY::OUTPUT"), vpflow.flag("SOLIQUIDY::X3D")).toString() << endl;
    } else {
      const xstructure input_structure(input, IOAFLOW_AUTO);
      const std::string result_folder = vpflow.getattachedscheme("SOLIQUIDY::OUTPUT");
      const aurostd::JSON::object result = Calculate(input_structure, result_folder, vpflow.flag("SOLIQUIDY::X3D"));
      const bool write_output = !result_folder.empty();
      const std::filesystem::path out_path(result_folder);
      if (write_output) {
        std::error_code err;
        if (!std::filesystem::create_directories(out_path, err)) {
          if (err.value() != 0) { // exist is ok
            const string message = "Can not create result folder " + result_folder + ". Error: " + err.message();
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
          }
        }
      }
      if (write_output) {
        result.saveFile((out_path / "soliquidy.json"));
      }
      cout << result.toString() << endl;
    }
  }

  /// @brief Load structures based on AUID and calculate their soliquidy value
  /// @param auid either one or vector of AUIDs
  /// @param result_folder folder to save the summary soliquidy.json and x3D content if requested (will be created if not exist)
  /// @param create_x3d flag to switch-on x3D creation
  /// @return results as JSON object
  template <class utype> aurostd::JSON::object CalculateAUID(const utype& auid, const std::string& result_folder, bool create_x3d) {
    const aflowlib::_aflowlib_entry entry;
    const xstructure str;
    aflowlib::EntryLoader el;
    el.loadAUID(auid);
    return CalculateEntryLoader(el.m_entries_flat, result_folder, create_x3d);
  }

  template aurostd::JSON::object CalculateAUID(const std::string&, const std::string&, bool);
  template aurostd::JSON::object CalculateAUID(const vector<std::string>&, const std::string&, bool);

  /// @brief Calculate Soliquidy for a list of structures, given by their file paths
  /// @param worklist vector of structure file paths
  /// @param result_folder folder to save the summary soliquidy.json and x3D content if requested (will be created if not exist)
  /// @param create_x3d flag to switch-on x3D creation
  /// @return results as JSON object
  aurostd::JSON::object CalculateWorklist(const std::vector<std::string>& worklist, const std::string& result_folder, bool create_x3d) {
    aurostd::JSON::object save_json(aurostd::JSON::object_types::DICTIONARY);
    const bool write_output = !result_folder.empty();
    const std::filesystem::path out_path(result_folder);
    if (write_output) {
      std::error_code err;
      if (!std::filesystem::create_directories(out_path, err)) {
        if (err.value() != 0) { // exist is ok
          const string message = "Can not create result folder " + result_folder + ". Error: " + err.message();
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
        }
      }
    }

    for (const auto& structure_file_raw : worklist) {
      const std::filesystem::path structure_file(structure_file_raw);
      const xstructure work_structure(structure_file_raw, IOAFLOW_AUTO);

      const std::filesystem::path single_out_folder = out_path / structure_file.stem();
      const aurostd::JSON::object result = Calculate(work_structure, single_out_folder, create_x3d);
      result["original_path"] = absolute(structure_file).string();
      save_json[structure_file.stem()] = result;
    }

    if (write_output) {
      save_json.saveFile((out_path / "soliquidy.json"));
    }
    return save_json;
  }

  /// @brief Calculate Soliquidy for all results in an EntryLoader view
  /// @param aflux_result flat view of an EntryLoader result
  /// @param result_folder folder to save the summary soliquidy.json and x3D content if requested (will be created if not exist)
  /// @param create_x3d flag to switch-on x3D creation
  /// @return results as JSON object
  aurostd::JSON::object CalculateEntryLoader(const std::shared_ptr<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>& aflux_result, const std::string& result_folder, bool create_x3d) {
    aurostd::JSON::object save_json(aurostd::JSON::object_types::DICTIONARY);
    const bool write_output = !result_folder.empty();
    const std::filesystem::path out_path(result_folder);
    if (write_output) {
      std::error_code err;
      if (!std::filesystem::create_directories(out_path, err)) {
        if (err.value() != 0) // exist is okay
        {
          const string message = "Can not create result folder " + result_folder + ". Error: " + err.message();
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
        }
      }
    }

    xstructure work_structure;
    for (const auto& entry : *aflux_result) {
      work_structure.clear();
      pflow::loadXstructureLibEntry(*entry, work_structure);
      const std::filesystem::path single_out_folder = out_path / (entry->auid).substr(6);
      const aurostd::JSON::object result = Calculate(work_structure, single_out_folder, create_x3d);
      save_json[entry->auid] = result;
    }

    if (write_output) {
      save_json.saveFile((out_path / "soliquidy.json"));
    }
    return save_json;
  }

  /// Calculate Soliquidy for a single structure
  /// @param loaded_structure xstructure to analyze
  /// @param out_folder folder to save x3D content if requested (will be created if not exist)
  /// @param create_x3d flag to switch-on x3D creation
  /// @return result as JSON object
  aurostd::JSON::object Calculate(const xstructure& loaded_structure, const std::filesystem::path& out_folder, const bool create_x3d) {
    const bool debug = (false || XHOST.DEBUG);
    const bool write_output = !out_folder.empty();

    // create a sub folder for the x3D representations
    if (write_output && create_x3d) {
      std::error_code err;
      if (!std::filesystem::create_directories(out_folder, err)) {
        if (err.value() != 0) { // exist is ok
          const string message = "Can not create result folder " + out_folder.string() + ". Error: " + err.message();
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
        }
      }
    }

    // the x3D object needs to be created outside the if-clause, or it will be out of scope when using it later
    aurostd::x3DWriter w;

    // copy the loaded xstructure and formate it for the use with aflow_voro
    xstructure work_structure = loaded_structure;
    work_structure.ReScale(1.0);

    // ensure that the lattice is defined in the aflow standard notation
    // note: similar to the standard_form() function in ASE

    xmatrix<double> Q;
    xmatrix<double> R;
    aurostd::QRDecomposition_HouseHolder_TB(aurostd::trasp(work_structure.lattice), Q, R);
    traspInPlace(Q);
    traspInPlace(R);

    // correct signs
    for (int row_idx = R.lrows; row_idx <= R.urows; row_idx++) {
      if (R[row_idx][row_idx] < 0) {
        for (int column_idx = R.lcols; column_idx <= R.ucols; column_idx++) {
          Q[row_idx][column_idx] *= (-1.0);
        }
      }
    }

    // apply the transformation
    work_structure.TransformStructure(
        {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
    },
        Q);
    work_structure.BringInCell();

    // define work variables and constants
    uint atom_count = 0;
    int p_id = 0;
    double smallest_cell_threshold = 0.0;
    double smallest_cell_area = 0.0;
    double cell_volume = 0.0;
    double target_volume = 0.0;
    double gradiant_norm = 0.0;
    double gradiant_norm_sub_step = 0.0;
    double alpha = 1.0;
    double split_r = 0.0;
    double stdev = 0.0;
    double stdev_compare = 0.0;
    double perfect_weight = 0.0;
    const double tikhonov_lambda = 1E-6;
    const double stdev_cutoff = 1E-6;
    bool KMT_1;
    bool KMT_2;
    std::function<double(xvector<double>)> density;
    voro::voronoicell_neighbor c;

    // prepare the data collection JSON object
    aurostd::JSON::object result(aurostd::JSON::object_types::DICTIONARY);
    result["species"] = aurostd::JSON::object(aurostd::JSON::object_types::DICTIONARY);
    result["value"] = 0.0;

    // loop over the species to calculate the induvidual soliquidy value for each
    for (const std::string& species : work_structure.species) {
      result["species"][species] = aurostd::JSON::object(aurostd::JSON::object_types::DICTIONARY);
      std::vector<xvector<double>> vertices;

      // collect the atom positions of the given species
      for (const _atom& atom : work_structure.atoms) {
        if (atom.cleanname == species) {
          vertices.push_back({atom.cpos[1], atom.cpos[2], atom.cpos[3]});
        }
      }
      atom_count = vertices.size();
      if (debug) {
        cout << species << " | count: " << atom_count << endl;
      }

      // define the storage for the fitting data
      const xvector<double> volumes(atom_count - 1, 0); // start x vector from 0 to align it with the cell id p_id
      xvector<double> gradiant(atom_count - 1, 0);
      xvector<double> weights(atom_count - 1, 0);
      xvector<double> weights_test(atom_count - 1, 0);
      xvector<double> p(atom_count - 1, 0);
      xmatrix<double> laplacian(atom_count - 1, atom_count - 1, 0, 0);
      xmatrix<double> laplacian_inverse(atom_count - 1, atom_count - 1, 0, 0);
      xmatrix<double> I(atom_count - 1, atom_count - 1, 0, 0);
      const xvector<double> soliquidy_cell_results(atom_count, 1);
      vector<double> soliquidy_cell_values;
      vector<xvector<double>> normals;
      vector<xvector<double>> cell_points;
      vector<vector<uint>> cell_facets;
      vector<xvector<double>> plot_points;
      vector<vector<uint>> plot_facets;
      double x;
      double y;
      double z;
      I = aurostd::identity(I);
      weights += 1.0; // set initial weights to 1 (xvectors are 0 initialized)

      { // section::original | calculate the unaltered voronoi cells
        // voro::container_periodic_poly is kept in its own environment as container.clear() does not work reliable
        voro::container_periodic_poly container(work_structure.lattice[1][1], work_structure.lattice[2][1], work_structure.lattice[2][2], work_structure.lattice[3][1], work_structure.lattice[3][2],
                                                work_structure.lattice[3][3], 6, 6, 6, 8);
        for (size_t i_vert = 0; i_vert < atom_count; i_vert++) {
          container.put(i_vert, vertices[i_vert][1], vertices[i_vert][2], vertices[i_vert][3], weights[i_vert]);
        }
        // loop over all cells
        if (voro::c_loop_all_periodic cl(container); cl.start()) {
          do {
            if (container.compute_cell(c, cl)) {
              p_id = cl.pid();
              volumes[p_id] = c.volume();
            }
          } while (cl.inc());
        }
      } // section::original

      // collect ststistic from the initial run
      cell_volume = aurostd::sum(volumes);
      target_volume = cell_volume / atom_count;
      if (debug) {
        cout << "    target_volume: " << target_volume << endl;
      }
      stdev = 0.0;
      for (int idx = volumes.lrows; idx <= volumes.urows; idx++) {
        stdev += std::pow(volumes[idx] - target_volume, 2);
      }
      stdev = sqrt(stdev) / target_volume;
      // check if the original voronoi cells are already optimal (for example in the case of a single atom)
      if (stdev < stdev_cutoff) {
        if (debug) {
          cout << "   No optimization needed!" << endl;
        }
      } else {
        // Smallest cell threshold for KMT criterion #1 and #2
        smallest_cell_area = aurostd::min(volumes);
        smallest_cell_threshold = 0.5 * std::min(smallest_cell_area, cell_volume / atom_count);

        // Newton steps with adaptive step length
        uint iteration_count = 0;
        do {
          laplacian.clear();
          { // section::step | take a Newton step
            // voro::container_periodic_poly is kept in its own environment as container.clear() does not work reliable
            voro::container_periodic_poly container(work_structure.lattice[1][1], work_structure.lattice[2][1], work_structure.lattice[2][2], work_structure.lattice[3][1], work_structure.lattice[3][2],
                                                    work_structure.lattice[3][3], 6, 6, 6, 8);
            for (size_t i_vert = 0; i_vert < atom_count; i_vert++) {
              container.put(i_vert, vertices[i_vert][1], vertices[i_vert][2], vertices[i_vert][3], weights[i_vert]);
            }
            if (voro::c_loop_all_periodic cl(container); cl.start()) {
              do {
                if (container.compute_cell(c, cl)) {
                  // fill Laplacian matrix
                  p_id = cl.pid();
                  c.fill_laplacian(p_id, vertices, work_structure.lattice, laplacian);
                  volumes[p_id] = c.volume();
                }
              } while (cl.inc());
            }
          } // section::step

          // find step length (alpha)
          alpha = 1.0;
          gradiant = -(volumes - target_volume);
          gradiant_norm = aurostd::modulus(gradiant);
          iteration_count++;
          laplacian = laplacian + I * tikhonov_lambda;

          laplacian_inverse = aurostd::inverseByQR(laplacian);
          laplacian_inverse.shift(0, 0);
          p = (laplacian_inverse * gradiant);

          if (debug) {
            cout << "  STEP " << iteration_count << " START" << endl;
          }
          // Start the substeps
          for (uint sub_step = 1; sub_step <= 6; sub_step++) {
            weights_test = weights + alpha * p;
            double offset_factor = 0.0;
            for (size_t i_vert = 0; i_vert < atom_count; i_vert++) {
              if (weights_test[i_vert] < 0.05) {
                offset_factor = max((0.05 - weights_test[i_vert]), offset_factor);
              }
            } // this ensures that weights never go below zero (0.05)

            if (offset_factor > 0) {
              if (debug) {
                cout << "      offset_factor: " << offset_factor << endl;
              }
            }
            weights_test += offset_factor;

            { // section::sub_step
              // voro::container_periodic_poly is kept in its own environment as container.clear() does not work reliable
              voro::container_periodic_poly container(work_structure.lattice[1][1], work_structure.lattice[2][1], work_structure.lattice[2][2], work_structure.lattice[3][1], work_structure.lattice[3][2],
                                                      work_structure.lattice[3][3], 6, 6, 6, 8);
              // set test weights
              for (size_t i_vert = 0; i_vert < atom_count; i_vert++) {
                container.put(i_vert, vertices[i_vert][1], vertices[i_vert][2], vertices[i_vert][3], weights_test[i_vert]);
              }
              // calculate system
              if (voro::c_loop_all_periodic cl(container); cl.start()) {
                do {
                  if (container.compute_cell(c, cl)) {
                    p_id = cl.pid();
                    volumes[p_id] = c.volume();
                  }
                } while (cl.inc());
              }
            } // section::sub_step

            // collect sub set data
            smallest_cell_area = aurostd::min(volumes);
            gradiant = -(volumes - target_volume);
            gradiant_norm_sub_step = aurostd::modulus(gradiant);
            KMT_1 = smallest_cell_area > smallest_cell_threshold;
            KMT_2 = gradiant_norm_sub_step <= ((1.0 - 0.5 * alpha) * gradiant_norm);
            stdev_compare = 0.0;
            for (int idx = volumes.lrows; idx <= volumes.urows; idx++) {
              stdev_compare += std::pow(volumes[idx] - target_volume, 2);
            }
            stdev_compare = sqrt(stdev_compare) / target_volume;

            if (KMT_1 and KMT_2) {
              break;
            }
            if (debug) {
              cout << "      KMT_1 (cell area): " << smallest_cell_area << " > " << smallest_cell_threshold << endl;
            }
            if (debug) {
              cout << "      KMT_2 (gradient): " << gradiant_norm_sub_step << " <= " << (1.0 - 0.5 * alpha) * gradiant_norm << endl;
            }
            alpha = alpha / 2.0;
          }

          weights = weights_test;
          stdev = stdev_compare;
          if (debug) {
            cout << "  STEP " << iteration_count << " DONE - stdev: " << stdev << endl;
            cout << "    final volumes: " << volumes << endl;
            cout << "    final weights: " << weights << endl;
            cout << "    stdev: " << stdev << endl;
            cout << "  --------------------------------------------------------------------" << endl;
          }
        } while (stdev > stdev_cutoff and iteration_count < 100);
      }

      // prepare result
      if (write_output && create_x3d) {
        w.clear();
        w.scene_center = (work_structure.lattice.getcol(1) + work_structure.lattice.getcol(2) + work_structure.lattice.getcol(3)) / 3.0;
        w.addLatticeBox(work_structure.lattice, 0.1);
      }
      { // section::final
        // voro::container_periodic_poly is kept in its own environment as container.clear() does not work reliable
        voro::container_periodic_poly container(work_structure.lattice[1][1], work_structure.lattice[2][1], work_structure.lattice[2][2], work_structure.lattice[3][1], work_structure.lattice[3][2],
                                                work_structure.lattice[3][3], 6, 6, 6, 8);
        // set final weights
        for (size_t i_vert = 0; i_vert < atom_count; i_vert++) {
          container.put(i_vert, vertices[i_vert][1], vertices[i_vert][2], vertices[i_vert][3], weights[i_vert]);
        }

        split_r = std::pow((3 * cell_volume) / (4 * PI * atom_count), 1.0 / 3.0);
        perfect_weight = PI * std::pow(split_r, 4);
        if (write_output && create_x3d) {
          plot_facets.clear();
          plot_points.clear();
        }
        if (voro::c_loop_all_periodic cl(container); cl.start()) {
          do {
            if (container.compute_cell(c, cl)) {
              p_id = cl.pid();
              cl.pos(x, y, z);
              cell_points.clear();
              cell_facets.clear();
              c.get_vertex_facets(x, y, z, cell_points, cell_facets);
              density = [x, y, z](const xvector<double>& point) -> double { return aurostd::distance(point, xvector<double>{x, y, z}); };
              const double tmp = aurostd::convex_volume_integral(cell_points, cell_facets, density);
              soliquidy_cell_results[p_id + 1] = std::pow((tmp - perfect_weight), 2);
              soliquidy_cell_values.emplace_back(tmp);
              if (write_output && create_x3d) {
                normals.clear();
                sortPolyhedronFacets(cell_points, cell_facets, normals);
                w.addSphere({x, y, z}, 0.3, "bm_blue");
                w.addSpheres(cell_points, 0.05);
                w.joinFacets(plot_points, plot_facets, cell_points, cell_facets);
              }
            }
          } while (cl.inc());

          result["species"][species]["stdev_species"] = stdev;
          result["species"][species]["values_cell"] = soliquidy_cell_values;
          result["species"][species]["perfect_volume"] = target_volume;
          result["species"][species]["perfect_weight"] = perfect_weight;
          result["species"][species]["value"] = 100 * std::sqrt(aurostd::mean(soliquidy_cell_results)) / perfect_weight;
          result["value"] = static_cast<double>(result["value"]) + (static_cast<double>(result["species"][species]["value"]) * atom_count);
        }
      } // section::final
      if (write_output && create_x3d) {
        w.addConvexFacets(plot_points, plot_facets);
        aurostd::string2file(w.toHTML(), (out_folder / (species + ".html")));
        aurostd::string2file(w.toTachyon(), (out_folder / (species + ".ty")));
        w.ani_type = aurostd::x3DWriter::animation_format::MP4;
        w.animate(5.0, out_folder / ("animation_" + species), 30, false);
      }
    }
    result["value"] = static_cast<double>(result["value"]) / work_structure.atoms.size();

    if (debug) {
      cout << "    soliquidy: " << result["value"].toString() << endl;
      for (const auto& [species_name, values] : static_cast<aurostd::JSON::Dictionary>(result["species"])) {
        cout << "    " << species_name << ": " << values["value"].toString() << endl;
      }
    }
    return result;
  }
} // namespace soliquidy

// ***************************************************************************
// Implementing JSON SERIALIZABLE HERE
// ***************************************************************************

aurostd::JSON::object xstructure::serialize() const {
  return aurostd::JSON::object({AST_JSON_GETTER(JSON_xstructure_MEMBERS)});
}

xstructure xstructure::deserialize(const aurostd::JSON::object& jo) {
  AST_JSON_SETTER(JSON_xstructure_MEMBERS)
  return *this;
}

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
