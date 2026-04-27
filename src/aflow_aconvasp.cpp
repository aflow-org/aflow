// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo and Dane Morgan

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <deque>
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
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_pflow.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

using std::deque;
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

using aurostd::xmatrix;
using aurostd::xvector;

// ***************************************************************************
// Simple interfaces
// ***************************************************************************
int SignNoZero(const double& x) {
  return (int) aurostd::signnozero(x);
}
int Nint(const double& x) {
  return (int) aurostd::nint(x);
}
int Sign(const double& x) {
  return (int) aurostd::sign(x);
}

// ***************************************************************************
// Function PutInCompact
// ***************************************************************************
// Make a structure where all the atoms are all the
// atoms are mapped through the unit and neighbors cells
// to minimixe the shortest possible bond with an adjacent atom. (SC 6 Aug 04)
xstructure PutInCompact(const xstructure& a) {
  xstructure sstr = a;
  sstr.BringInCompact();
  return sstr;
}

// ***************************************************************************
// Function PutInCell - same as BringInCell
// ***************************************************************************
// Make a structure with all atoms mapped to their images within the
// unit cell.
xstructure PutInCell(const xstructure& a) {
  xstructure sstr = a;
  sstr.BringInCell();
  return sstr;
}

// ***************************************************************************
// Function WignerSeitz
// ***************************************************************************
// Make a structure where all the atoms are
// mapped to their images in the Wigner-Seitz cell.(SC 10Jan04)
xstructure WignerSeitz(const xstructure& a) {
  xstructure sstr = a;
  // zoology of vectors
  xvector<double> rat(3);
  xvector<double> xoo(3);
  xvector<double> yoo(3);
  xvector<double> zoo(3);
  xvector<double> xyo(3);
  xvector<double> xzo(3);
  xvector<double> yzo(3);
  xvector<double> xyz(3);
  double axoo;
  double ayoo;
  double azoo;
  double axyo;
  double axzo;
  double ayzo;
  double axyz;
  xoo = sstr.lattice(1);
  axoo = modulus(xoo);
  yoo = sstr.lattice(2);
  ayoo = modulus(yoo);
  zoo = sstr.lattice(3);
  azoo = modulus(zoo);
  xyo = sstr.lattice(1) + sstr.lattice(2);
  axyo = modulus(xyo);
  xzo = sstr.lattice(1) + sstr.lattice(3);
  axzo = modulus(xzo);
  yzo = sstr.lattice(2) + sstr.lattice(3);
  ayzo = modulus(yzo);
  xyz = sstr.lattice(1) + sstr.lattice(2) + sstr.lattice(3);
  axyz = modulus(xyz);
  double projxoo;
  double projyoo;
  double projzoo;
  double projxyo;
  double projxzo;
  double projyzo;
  double projxyz;

  for (size_t iat = 0; iat < sstr.atoms.size(); iat++) {
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          rat = sstr.atoms.at(iat).cpos + ((double) i) * sstr.lattice(1) + ((double) j) * sstr.lattice(2) + ((double) k) * sstr.lattice(3);
          projxoo = scalar_product(rat, xoo) / axoo / axoo;
          projyoo = scalar_product(rat, yoo) / ayoo / ayoo;
          projzoo = scalar_product(rat, zoo) / azoo / azoo;
          projxyo = scalar_product(rat, xyo) / axyo / axyo;
          projxzo = scalar_product(rat, xzo) / axzo / axzo;
          projyzo = scalar_product(rat, yzo) / ayzo / ayzo;
          projxyz = scalar_product(rat, xyz) / axyz / axyz;
          if ((projxoo > -0.5 && projxoo <= 0.5) && (projyoo > -0.5 && projyoo <= 0.5) && (projzoo > -0.5 && projzoo <= 0.5) && (projxyo > -0.5 && projxyo <= 0.5) && (projxzo > -0.5 && projxzo <= 0.5) &&
              (projyzo > -0.5 && projyzo <= 0.5) && (projxyz > -0.5 && projxyz <= 0.5)) {
            sstr.atoms.at(iat).cpos(1) = rat(1);
            sstr.atoms.at(iat).cpos(2) = rat(2);
            sstr.atoms.at(iat).cpos(3) = rat(3);
            i = 10;
            j = 10;
            k = 10;
          }
        }
      }
    }
    sstr.atoms.at(iat).fpos = C2F(sstr.lattice, sstr.atoms.at(iat).cpos);
  }
  return sstr;
}

// **************************************************************************
// GetMom1
// **************************************************************************
// this function returns moment_1 position of the atoms in cartesians
xvector<double> GetMom1(const xstructure& a) {
  // Get's first moment in cartesian coordinates.
  xvector<double> mom1(3);
  for (size_t iat = 0; iat < a.atoms.size(); iat++) {
    mom1 = mom1 + a.atoms[iat].cpos;
  }
  if (!a.atoms.empty()) {
    mom1 = ((double) a.scale / ((double) a.atoms.size())) * mom1;
  } else {
    clear(mom1);
  }
  return mom1;
}

// **************************************************************************
// SetMom1
// **************************************************************************
// this function sets moment_1 position of atoms
xstructure SetMom1(const xstructure& a, const xvector<double>& mom1_new) {
  xstructure sstr = a;
  // Get old first moment in cart. coords.
  xvector<double> mom1_old(3);
  mom1_old = GetMom1(sstr);
  // Get change in first moment in cart. coords.
  xvector<double> dmom1(3);
  dmom1 = mom1_new - mom1_old;
  // Shift scaled cart. positions by dmom1.
  for (size_t iat = 0; iat < sstr.atoms.size(); iat++) {
    sstr.atoms[iat].cpos = sstr.scale * sstr.atoms[iat].cpos;
    sstr.atoms[iat].cpos = sstr.atoms[iat].cpos + dmom1;
    sstr.atoms[iat].cpos = sstr.atoms[iat].cpos / sstr.scale;
    sstr.atoms[iat].fpos = C2F(sstr.lattice, sstr.atoms[iat].cpos);
  }
  return sstr;
}

// **************************************************************************
// Function Function GetCDispFromOrigin
// **************************************************************************
// This function returns displacement from origin
xvector<double> GetCDispFromOrigin(const _atom& atom) {
  xvector<double> diff(3);
  diff = atom.cpos - atom.corigin;
  return diff;
}

// **************************************************************************
// Function GetDistFromOrigin
// **************************************************************************
// This function returns distance from origin
double GetDistFromOrigin(const _atom& atom) {
  return modulus(GetCDispFromOrigin(atom));
}

// **************************************************************************
// Function ConvertAtomToLat
// **************************************************************************
// This function makes sure that an atoms fractional coords
// and unit cell parameters are consistent with an input lattice.
_atom ConvertAtomToLat(const _atom& in_at, const xmatrix<double>& lattice) {
  xvector<double> cpos(3);
  xvector<double> p_cell0(3);
  xvector<int> ijk;
  _atom out_at = in_at;
  cpos = in_at.cpos;
  GetUnitCellRep(cpos, p_cell0, ijk, lattice, true);
  out_at.fpos = C2F(lattice, cpos);
  out_at.ijk = ijk;
  return out_at;
}

// **************************************************************************
// Function GetXrayScattFactor
// **************************************************************************
double GetXrayScattFactor(const string& _name, double lambda, bool clean) {
  string name = _name;  // CO20190322
  if (clean) {
    name = aurostd::VASP_PseudoPotential_CleanName(name);
  } // CO20190322
  if (lambda) {
    ;
  } // phony just to keep lambda busy
  // Does not use lambda for now.
  double scatt_fact = 0.0;
  for (size_t iat = 0; iat < vatom_name.size(); iat++) {
    if (name == vatom_name[iat] || name == vatom_symbol[iat]) {
      scatt_fact = vatom_xray_scatt[iat];
    }
  }
  return scatt_fact;
}

// ***************************************************************************
// Function GetVol
// ***************************************************************************
// This function returns the volume of the cell from the lattice parameters
namespace pflow {
  double GetVol(const xmatrix<double>& lat) {
    return aurostd::abs(pflow::GetSignedVol(lat));
  }
  double GetVol(const aurostd::matrix<double>& lat) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    return aurostd::abs(pflow::GetSignedVol(lat));
  }
} // namespace pflow

// ***************************************************************************
// Function GetSignedVol
// ***************************************************************************
// This function returns the volume of the cell from the lattice parameters
namespace pflow {
  double GetSignedVol(const xmatrix<double>& lat) {
    double vol;
    xvector<double> u(3);
    u = aurostd::vector_product(lat(2), lat(3));
    vol = aurostd::scalar_product(lat(1), u);
    return vol;
  }
  double GetSignedVol(const aurostd::matrix<double>& lat) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    return pflow::GetSignedVol(aurostd::matrix2xmatrix(lat)); // CO20200404 pflow::matrix()->aurostd::matrix()
  }
} // namespace pflow

// ***************************************************************************
// Function RecipLat xmatrix<double>
// ***************************************************************************
namespace pflow {
  xmatrix<double> RecipLat(const xmatrix<double>& lat) {
    xmatrix<double> rlat(3, 3);
    xvector<double> rlat1(3);
    xvector<double> rlat2(3);
    xvector<double> rlat3(3);
    const double vol = pflow::GetSignedVol(lat);
    rlat1 = (2.0 * PI / vol) * aurostd::vector_product(lat(2), lat(3));
    rlat2 = (2.0 * PI / vol) * aurostd::vector_product(lat(3), lat(1));
    rlat3 = (2.0 * PI / vol) * aurostd::vector_product(lat(1), lat(2));
    for (int i = 1; i <= 3; i++) {
      rlat(1, i) = rlat1(i);
    }
    for (int i = 1; i <= 3; i++) {
      rlat(2, i) = rlat2(i);
    }
    for (int i = 1; i <= 3; i++) {
      rlat(3, i) = rlat3(i);
    }
    return rlat;
  }
} // namespace pflow

// ***************************************************************************
// Function RecipLat aurostd::matrix<double>
// ***************************************************************************
namespace pflow {
  aurostd::matrix<double> RecipLat(const aurostd::matrix<double>& lat) { // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rlat(3, 3);  // CO20200404 pflow::matrix()->aurostd::matrix()
    const double vol = pflow::GetSignedVol(lat);
    rlat[0] = pflow::SVprod(2.0 * PI / vol, pflow::VVcross(lat[1], lat[2]));
    rlat[1] = pflow::SVprod(2.0 * PI / vol, pflow::VVcross(lat[2], lat[0]));
    rlat[2] = pflow::SVprod(2.0 * PI / vol, pflow::VVcross(lat[0], lat[1]));
    return rlat;
  }
} // namespace pflow

// ***************************************************************************
// Function SetCpos
// ***************************************************************************
namespace pflow {
  _atom SetCpos(const _atom& a, const vector<double>& in_cpos) {
    assert(in_cpos.size() == 3);
    _atom b;
    b = a;
    b.cpos(1) = in_cpos.at(0);
    b.cpos(2) = in_cpos.at(1);
    b.cpos(3) = in_cpos.at(2);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function SetFpos
// ***************************************************************************
namespace pflow {
  _atom SetFpos(const _atom& a, const vector<double>& in_fpos) {
    assert(in_fpos.size() == 3);
    _atom b;
    b = a;
    b.fpos(1) = in_fpos.at(0);
    b.fpos(2) = in_fpos.at(1);
    b.fpos(3) = in_fpos.at(2);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function vecF2C
// ***************************************************************************
// This function converts a vector in direct to cartesian
namespace pflow {
  vector<double> vecF2C(const aurostd::matrix<double>& lat, const vector<double>& vf) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> vc(3, 0.0);
    for (int ic = 0; ic < 3; ic++) {
      for (int jc = 0; jc < 3; jc++) {
        vc[ic] = vc[ic] + vf[jc] * lat[jc][ic];
      }
    }
    return vc;
  }
} // namespace pflow

// ***************************************************************************
// Function vecC2F
// ***************************************************************************
// This function converts a vector in direct to cartesian
namespace pflow {
  vector<double> vecC2F(const aurostd::matrix<double>& lat, const vector<double>& vc) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> vf(3, 0.0);
    vf = aurostd::xvector2vector(C2F(matrix2xmatrix(lat), aurostd::vector2xvector(vc)));
    return vf;
  }
} // namespace pflow

// ***************************************************************************
// Function SetName
// ***************************************************************************
namespace pflow {
  _atom SetName(const _atom& a, const string& in_name) {
    _atom b;
    b = a;
    b.name = in_name;
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function SetType
// ***************************************************************************
namespace pflow {
  _atom SetType(const _atom& a, const int in_type) {
    _atom b;
    b = a;
    b.type = in_type;    // CONVASP_MODE
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function SetNum
// ***************************************************************************
namespace pflow {
  _atom SetNum(const _atom& a, const int in_num) {
    _atom b;
    b = a;
    b.basis = in_num;  //[CO20200130 - number->basis]b.number=in_num;
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function GetFpos
// ***************************************************************************
namespace pflow {
  aurostd::matrix<double> GetFpos(const xstructure& str) { // CO20200404 pflow::matrix()->aurostd::matrix()
    const int num_atoms = str.atoms.size();
    aurostd::matrix<double> fpos(num_atoms, 3);
    pflow::VVset(fpos, 0.0);  // CO20200404 pflow::matrix()->aurostd::matrix()
    for (int i = 0; i < num_atoms; i++) {
      for (int j = 0; j < 3; j++) {
        fpos[i][j] = str.atoms.at(i).fpos(j + 1);
      }
    }
    return fpos;
  }
} // namespace pflow

// ***************************************************************************
// Function GetCpos
// ***************************************************************************
namespace pflow {
  aurostd::matrix<double> GetCpos(const xstructure& str) { // CO20200404 pflow::matrix()->aurostd::matrix()
    const int num_atoms = str.atoms.size();
    aurostd::matrix<double> cpos(num_atoms, 3);
    pflow::VVset(cpos, 0.0);  // CO20200404 pflow::matrix()->aurostd::matrix()
    for (int i = 0; i < num_atoms; i++) {
      for (int j = 0; j < 3; j++) {
        cpos[i][j] = str.atoms.at(i).cpos(j + 1);
      }
    }
    return cpos;
  }
} // namespace pflow

// ***************************************************************************
// Function SetLat
// ***************************************************************************
xstructure SetLat(const xstructure& a, const xmatrix<double>& in_lat) {
  xstructure b(a);
  b.lattice = in_lat;
  return b;
}
namespace pflow {
  xstructure SetLat(const xstructure& a, const aurostd::matrix<double>& in_lat) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    xstructure b(a);
    b.lattice = aurostd::matrix2xmatrix(in_lat); // CO20200404 pflow::matrix()->aurostd::matrix()
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function GetLat
// ***************************************************************************
xmatrix<double> GetLat(const xstructure& a) {
  return a.lattice;
}
namespace pflow {
  aurostd::matrix<double> GetLat(const xstructure& a) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    return aurostd::xmatrix2matrix(a.lattice); // CO20200404 pflow::matrix()->aurostd::matrix()
  }
} // namespace pflow

// ***************************************************************************
// Function GetScale
// ***************************************************************************
namespace pflow {
  double GetScale(const xstructure& a) {
    return a.scale;
  }
} // namespace pflow

// ***************************************************************************
// Function GetScaledLat
// ***************************************************************************
namespace pflow {
  aurostd::matrix<double> GetScaledLat(const xstructure& a) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    return aurostd::xmatrix2matrix(a.scale * a.lattice); // CO20200404 pflow::matrix()->aurostd::matrix()
  }
} // namespace pflow

// ***************************************************************************
// ***************************************************************************

// ***************************************************************************
// ***************************************************************************

// ***************************************************************************
// Function AddAllAtomPos
// ***************************************************************************
namespace pflow {
  xstructure AddAllAtomPos(const xstructure& a,
                           const aurostd::matrix<double>& in_pos, // CO20200404 pflow::matrix()->aurostd::matrix()
                           const int in_coord_flag) {
    assert(in_coord_flag == 0 || in_coord_flag == 1);
    xstructure b(a);
    b.atoms.clear();
    _atom atom;
    if (in_coord_flag == 0) {
      for (size_t i = 0; i < in_pos.size(); i++) {
        // atom=a.atoms.at(i);  // start from scratch
        for (int j = 1; j <= 3; j++) {
          atom.fpos(j) = in_pos[i][j - 1];
        }
        atom.cpos = F2C(b.lattice, atom.fpos);
        b.atoms.push_back(atom);
      }
    }
    if (in_coord_flag == 1) {
      for (size_t i = 0; i < in_pos.size(); i++) {
        // atom=a.atoms.at(i);  // start from scratch
        for (int j = 1; j <= 3; j++) {
          atom.cpos(j) = in_pos[i][j - 1];
        }
        atom.fpos = C2F(b.lattice, atom.cpos);
        b.atoms.push_back(atom);
      }
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function SetAllAtomPos
// ***************************************************************************
namespace pflow {
  xstructure SetAllAtomPos(const xstructure& a,
                           const aurostd::matrix<double>& in_pos, // CO20200404 pflow::matrix()->aurostd::matrix()
                           const int in_coord_flag) {
    assert(in_coord_flag == 0 || in_coord_flag == 1);
    xstructure b(a);
    b.atoms.clear();
    _atom atom;
    if (in_coord_flag == 0) {
      for (size_t i = 0; i < in_pos.size(); i++) {
        atom = a.atoms.at(i);
        for (int j = 1; j <= 3; j++) {
          atom.fpos(j) = in_pos[i][j - 1];
        }
        atom.cpos = F2C(b.lattice, atom.fpos);
        b.atoms.push_back(atom);
      }
    }
    if (in_coord_flag == 1) {
      for (size_t i = 0; i < in_pos.size(); i++) {
        atom = a.atoms.at(i);
        for (int j = 1; j <= 3; j++) {
          atom.cpos(j) = in_pos[i][j - 1];
        }
        atom.fpos = C2F(b.lattice, atom.cpos);
        b.atoms.push_back(atom);
      }
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function SetAllAtomNames
// ***************************************************************************
namespace pflow {
  xstructure SetAllAtomNames(const xstructure& a, const vector<string>& in) {
    xstructure b(a);
    if (in.size() == a.num_each_type.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "this routine must be fixed, it does not work here", _RUNTIME_ERROR_);

      for (size_t iat = 0; iat < b.num_each_type.size(); iat++) {
        b.atoms.at(iat).name = in.at(b.atoms.at(iat).type);       // CONVASP_MODE
        if (!b.atoms.at(iat).name.empty()) {
          b.atoms.at(iat).name_is_given = true;
          b.atoms.at(iat).CleanName();
          // DX20170921 - Need to keep spin info b.atoms.at(iat).CleanSpin();
        }
      }
      return b;
    }
    if (in.size() == a.atoms.size()) {
      for (size_t iat = 0; iat < b.atoms.size(); iat++) {
        b.atoms[iat].name = in[iat];
        b.atoms[iat].name_is_given = true;
        b.atoms[iat].CleanName();
        // DX20170921 - Need to keep spin info b.atoms.at(iat).CleanSpin();
      }
      return b;
    }
    stringstream message;
    message << "Must specify as many names as types/bases: in.size()=" << in.size();  //[CO20200130 - number->basis]
    message << "   =a.num_each_type.size()=" << a.num_each_type.size();
    message << "   =a.atoms.size()=" << a.atoms.size();
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
  }
} // namespace pflow

// ***************************************************************************
// Function SetNamesWereGiven
// ***************************************************************************
namespace pflow {
  xstructure SetNamesWereGiven(const xstructure& a, const vector<int>& in) {
    xstructure b(a);
    if (in.size() == a.num_each_type.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "this routine must be fixed... it does not work here", _RUNTIME_ERROR_);
      for (size_t iat = 0; iat < b.num_each_type.size(); iat++) {
        b.atoms.at(iat).name_is_given = in.at(b.atoms.at(iat).type);      // CONVASP_MODE
      }
      return b;
    }
    if (in.size() == a.atoms.size()) {
      for (size_t iat = 0; iat < b.atoms.size(); iat++) {
        b.atoms[iat].name_is_given = in[iat];
      }
      return b;
    }
    stringstream message;
    message << "Must specify as many names as types/bases: in.size()=" << in.size();  //[CO20200130 - number->basis]
    message << "   =a.num_each_type.size()=" << a.num_each_type.size();
    message << "   =a.atoms.size()=" << a.atoms.size();
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
  }
} // namespace pflow

// ***************************************************************************
// Function SetOrigin
// ***************************************************************************
namespace pflow {
  xstructure SetOrigin(const xstructure& a, const vector<double>& in_origin) {
    xstructure b(a);
    for (uint i = 1; i <= 3; i++) {
      b.origin[i] = in_origin.at(i - 1);
    }
    return b;
  }
  xstructure SetOrigin(const xstructure& a, const xvector<double>& in_origin) {
    xstructure b(a);
    for (uint i = 1; i <= 3; i++) {
      b.origin[i] = in_origin[i];
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// Function VVequal
// ***************************************************************************
// This function checks if two vectors are equal. double/int
namespace pflow {
  // vector
  bool VVequal(const vector<double>& a, const vector<double>& b) {
    const double tol = 1e-15;
    const int size = std::min((int) a.size(), (int) b.size());
    for (int i = 0; i < size; i++) {
      if (aurostd::abs(a[i] - b[i]) > tol) {
        return false;
      }
    }
    return true;
  }
  bool VVequal(const vector<int>& a, const vector<int>& b) {
    const double tol = 1e-15;
    const int size = std::min((int) a.size(), (int) b.size());
    for (int i = 0; i < size; i++) {
      if (aurostd::abs((double) (a[i] - b[i])) > tol) {
        return false;
      }
    }
    return true;
  }
  // deque
  bool VVequal(const deque<double>& a, const deque<double>& b) {
    const double tol = 1e-15;
    const int size = std::min((int) a.size(), (int) b.size());
    for (int i = 0; i < size; i++) {
      if (aurostd::abs(a[i] - b[i]) > tol) {
        return false;
      }
    }
    return true;
  }
  bool VVequal(const deque<int>& a, const deque<int>& b) {
    const double tol = 1e-15;
    const int size = std::min((int) a.size(), (int) b.size());
    for (int i = 0; i < size; i++) {
      if (aurostd::abs((double) (a[i] - b[i])) > tol) {
        return false;
      }
    }
    return true;
  }

} // namespace pflow
// ***************************************************************************
// Function SmoothFunc
// ***************************************************************************
// This function smooths an input functions.
// Uses gaussian smoothing for now.
// Assumes function argument is integer given by vector id.
// Sigma is in units of function arguments (i.e., vector id).
// Dane Morgan
namespace pflow {
  vector<double> SmoothFunc(const vector<double>& func, const double& sigma) {
    if (sigma < 1e-10) {
      return func;
    }
    const int nx = func.size();
    const int range = (int) (5.0 * sigma);
    // Get Gaussian weights
    vector<double> norm(nx, 0.0);
    aurostd::matrix<double> wt(nx, 2 * range + 1);
    pflow::VVset(wt, 0.0); // CO20200404 pflow::matrix()->aurostd::matrix()

    for (int ix = 0; ix < nx; ix++) {
      for (int i = -range; i <= range; i++) {
        if ((ix + i) >= 0 && (ix + i) < nx) {
          wt[ix][i + range] = Normal((double) (ix + i), double(ix), sigma);
          norm[ix] = norm[ix] + wt[ix][i + range];
        }
      }
    }
    // Normalize to one
    for (int ix = 0; ix < nx; ix++) {
      for (int i = -range; i <= range; i++) {
        wt[ix][i + range] = wt[ix][i + range] / norm[ix];
      }
    }

    vector<double> sfunc(nx, 0.0);
    // Average in weighted nearby bins.
    for (int ix = 0; ix < nx; ix++) {
      double sf = 0;
      for (int i = -range; i <= range; i++) {
        if ((ix + i) > 0 && (ix + i) < nx) {
          sf += wt[ix][i + range] * func[ix + i];
        }
      }
      sfunc[ix] = sf;
    }
    return sfunc;
  }
} // namespace pflow
// ***************************************************************************
//  Function Normal
// ***************************************************************************
// This function returns the value of a normal distribution.
double Normal(const double& x, const double& mu, const double& sigma) {
  const double tol = 1e-12;
  // If sigma=0 return delta function
  if (std::abs(sigma) < tol) {
    if (std::abs(x - mu) < tol) {
      return 1;
    } else {
      return 0;
    }
  }
  const double arg = -(x - mu) * (x - mu) / (2 * sigma * sigma);
  return exp(arg) / (sqrt(TWOPI) * sigma);
}

// ***************************************************************************
// Function Set
// ***************************************************************************
// This function sets a vector of vector with the proper values
namespace pflow {
  void VVset(aurostd::matrix<double>& mat, const double& value) { // CO20200404 pflow::matrix()->aurostd::matrix()
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[i].size(); j++) {
        mat[i][j] = (double) value;
      }
    }
  }
  void VVset(vector<vector<int>>& mat, const int& value) {
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[i].size(); j++) {
        mat[i][j] = (int) value;
      }
    }
  }
} // namespace pflow

// ***************************************************************************
// Function norm
// ***************************************************************************
// This function finds the norm of a vector
// Dane Morgan style
namespace pflow {
  //  template<class utype>
  double norm(const vector<double>& v) {
    double sum = 0;
    for (size_t i = 0; i < v.size(); i++) {
      sum = sum + v[i] * v[i];
    }
    //    sum=sqrt(aurostd::abs(sum));
    sum = sqrt(aurostd::abs(sum));
    return sum;
  }
} // namespace pflow

// ***************************************************************************
// Function getcos
// ***************************************************************************
// This function finds the cosine of the angle between
// two vectors.
// Dane Morgan style
namespace pflow {
  // template<class utype>
  double getcos(const vector<double>& a, const vector<double>& b) {
    double sum = 0.0;
    const uint size = std::min((uint) a.size(), (uint) b.size());
    for (uint i = 0; i < size; i++) {
      sum = sum + a[i] * b[i];
    }
    const double na = pflow::norm(a);
    const double nb = pflow::norm(b);
    assert(na > 0.0 && nb > 0.0);
    sum = sum / (na * nb);
    return sum;
  }
} // namespace pflow

// ***************************************************************************
// Function Getabc_angles
// ***************************************************************************
// This function returns a,b,c,alpha,beta,gamma for the cell
// given the lattice vectors.
// Dane Morgan style
namespace pflow {
  // template<class utype>
  vector<double> Getabc_angles(const aurostd::matrix<double>& lat) { // CO20200404 pflow::matrix()->aurostd::matrix()
    //    cerr << lat[0][0] << " " << lat[0][1] << " " << lat[0][2] << endl;
    //    cerr << lat[1][0] << " " << lat[1][1] << " " << lat[1][2] << endl;
    //   cerr << lat[2][0] << " " << lat[2][1] << " " << lat[2][2] << endl;
    vector<double> data;
    data.push_back(pflow::norm(lat[0])); // cerr << data.at(data.size()-1) << endl;
    data.push_back(pflow::norm(lat[1])); // cerr << data.at(data.size()-1) << endl;
    data.push_back(pflow::norm(lat[2])); // cerr << data.at(data.size()-1) << endl;
    data.push_back(acos(pflow::getcos(lat[1], lat[2]))); // cerr << data.at(data.size()-1) << endl;
    data.push_back(acos(pflow::getcos(lat[0], lat[2]))); // cerr << data.at(data.size()-1) << endl;
    data.push_back(acos(pflow::getcos(lat[0], lat[1]))); // cerr << data.at(data.size()-1) << endl;
    return data;
  }
} // namespace pflow

// namespace pflow {
//   void dont_run_this(void) {
//     {
//       vector<float> d;pflow::norm(d);pflow::getcos(d,d);
//       aurostd::matrix<float> D;pflow::Getabc_angles(D); //CO20200404 pflow::matrix()->aurostd::matrix()
//     }
//     {
//       vector<double> d;pflow::norm(d);pflow::getcos(d,d);
//       aurostd::matrix<double> D;pflow::Getabc_angles(D);  //CO20200404 pflow::matrix()->aurostd::matrix()
//     }
//     {
//       vector<long double> d;pflow::norm(d);pflow::getcos(d,d);
//       aurostd::matrix<long double> D;pflow::Getabc_angles(D); //CO20200404 pflow::matrix()->aurostd::matrix()
//     }
//   }
// }

// ***************************************************************************
// Function Sort_abc_angles
// ***************************************************************************
// This fucntions sorts lattice vectors to a>=b>=c and
// arranges angles equivalently.
// Dane Morgan style
namespace pflow {
  vector<double> Sort_abc_angles(const vector<double>& abc_angles) {
    vector<double> old = abc_angles;
    vector<double> newv = abc_angles;
    // if b>a then swap a/b: a=b,b=a,c=c,alpha=beta,beta=alpha,gamma=gamma.
    if (old[1] > old[0]) {
      newv[0] = old[1];
      newv[1] = old[0];
      newv[3] = old[4];
      newv[4] = old[3];
      old = newv;
    }
    // if c>a then swap a/c: a=c,b=b,c=a,alpha=gamma,beta=beta,gamma=alpha.
    if (old[2] > old[0]) {
      newv[0] = old[2];
      newv[2] = old[0];
      newv[3] = old[5];
      newv[5] = old[3];
      old = newv;
    }
    // if c>b then swap b/c: a=a,b=c,c=b,alpha=alpha,beta=gamma,gamma=beta.
    if (old[2] > old[1]) {
      newv[1] = old[2];
      newv[2] = old[1];
      newv[4] = old[5];
      newv[5] = old[4];
      old = newv;
    }
    return newv;
  }
} // namespace pflow

// ***************************************************************************
// Function Vout
// ***************************************************************************
// This function outputs a vector.
// Dane Morgan style
namespace pflow {
  void Vout(const vector<double>& a, ostream& out) {
    for (size_t i = 0; i < a.size(); i++) {
      out << a[i] << " ";
    }
    out << endl;
  }
  void Vout(const vector<int>& a, ostream& out) {
    for (size_t i = 0; i < a.size(); i++) {
      out << a[i] << " ";
    }
    out << endl;
  }
  void Vout(const vector<string>& a, ostream& out) {
    for (size_t i = 0; i < a.size(); i++) {
      out << a[i] << " ";
    }
    out << endl;
  }
} // namespace pflow

// ***************************************************************************
//  Function Mout
// ***************************************************************************
// This function outputs a matrix.
// Dane Morgan style
namespace pflow {
  void Mout(const aurostd::matrix<double>& m, ostream& out) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    for (size_t i = 0; i < m.size(); i++) {
      for (size_t j = 0; j < m[i].size(); j++) {
        out << m[i][j] << " ";
      }
      out << endl;
    }
  }
  void Mout(const vector<vector<double>>& m, ostream& out) {
    for (size_t i = 0; i < m.size(); i++) {
      for (size_t j = 0; j < m[i].size(); j++) {
        out << m[i][j] << " ";
      }
      out << endl;
    }
  }
} // namespace pflow

// **************************************************
//  Function SVprod
// **************************************************
// This function returns the product of a scalar and a vector.
namespace pflow {
  vector<double> SVprod(const double& s, const vector<double>& b) {
    const int size = b.size();
    vector<double> prod(size);
    for (int i = 0; i < size; i++) {
      prod[i] = s * b[i];
    }
    return prod;
  }
  vector<int> SVprod(const int& s, const vector<int>& b) {
    const int size = b.size();
    vector<int> prod(size);
    for (int i = 0; i < size; i++) {
      prod[i] = s * b[i];
    }
    return prod;
  }
} // namespace pflow

// **************************************************
//  Function VVsum
// **************************************************
// This function returns the vector sum c=a+b.
namespace pflow {
  vector<double> VVsum(const vector<double>& a, const vector<double>& b) {
    const int size = std::min(a.size(), b.size());
    vector<double> c(size);
    for (int i = 0; i < size; i++) {
      c[i] = a[i] + b[i];
    }
    return c;
  }
  vector<double> VVsum(const vector<double>& a, const vector<int>& b) {
    const int size = std::min(a.size(), b.size());
    vector<double> c(size);
    for (int i = 0; i < size; i++) {
      c[i] = a[i] + b[i];
    }
    return c;
  }
} // namespace pflow

// **************************************************
// Function VVdiff
// **************************************************
// This function returns the vector c=a-b.
namespace pflow {
  vector<double> VVdiff(const vector<double>& a, const vector<double>& b) {
    const int size = std::min(a.size(), b.size());
    vector<double> c(size);
    for (int i = 0; i < size; i++) {
      c[i] = a[i] - b[i];
    }
    return c;
  }
} // namespace pflow

// **************************************************
//  Function VVprod
// **************************************************
// This function returns the scalar c=a*b.
// Dane Morgan style
namespace pflow {
  double VVprod(const vector<double>& a, const vector<double>& b) {
    const int size = std::min(a.size(), b.size());
    double c = 0;
    for (int i = 0; i < size; i++) {
      c = c + a[i] * b[i];
    }
    return c;
  }
  double VVprod(const vector<double>& a, const vector<int>& b) {
    const int size = std::min(a.size(), b.size());
    double c = 0;
    for (int i = 0; i < size; i++) {
      c = c + a[i] * b[i];
    }
    return c;
  }
} // namespace pflow

// ***************************************************************************
// Function MMmult
// ***************************************************************************
// This function returns the product of two matrices.
// a=MxN, b=NxM.
// Dane Morgan style
namespace pflow {
  aurostd::matrix<double> MMmult(const aurostd::matrix<double>& a, const aurostd::matrix<double>& b) { // CO20200404 pflow::matrix()->aurostd::matrix()
    const uint M = std::min((uint) a.size(), (uint) b[0].size());
    const uint N = std::min((uint) a[0].size(), (uint) b.size());
    const vector<double> v(M, 0.0);
    aurostd::matrix<double> c(M, v);  // CO20200404 pflow::matrix()->aurostd::matrix()
    for (uint i = 0; i < M; i++) {
      for (uint j = 0; j < M; j++) {
        double sum = 0.0;
        for (uint k = 0; k < N; k++) {
          sum = sum + a[i][k] * b[k][j];
        }
        c[i][j] = sum;
      }
    }
    return c;
  }
} // namespace pflow

// ***************************************************************************
// Function MVmult
// ***************************************************************************
//  This function returns the product of a matric and a vector A*v.
//  A=MxN, v=Nx1.
// Dane Morgan style
namespace pflow {
  vector<double> MVmult(const aurostd::matrix<double>& A, const vector<double>& v) { // CO20200404 pflow::matrix()->aurostd::matrix()
    const uint M = A.size();
    const uint N = std::min((uint) A[0].size(), (uint) v.size());
    vector<double> u(M, 0.0);
    for (uint i = 0; i < M; i++) {
      double sum = 0;
      for (uint j = 0; j < N; j++) {
        sum = sum + A[i][j] * v[j];
      }
      u[i] = sum;
    }
    return u;
  }
} // namespace pflow

// ***************************************************************************
// Function VMmult
// ***************************************************************************
// This function returns the product of a row vector and a matrix v*A.
// A=MxN, v=1xM. Return vector of length N.
// Dane Morgan style
namespace pflow {
  vector<double> VMmult(const vector<double>& v, const aurostd::matrix<double>& A) { // CO20200404 pflow::matrix()->aurostd::matrix()
    const uint M = A.size();
    const uint N = std::min((uint) A[0].size(), (uint) v.size());
    vector<double> u(M, 0.0);
    for (uint i = 0; i < N; i++) {
      double sum = 0;
      for (uint j = 0; j < M; j++) {
        sum = sum + v[j] * A[j][i];
      }
      u[i] = sum;
    }
    return u;
  }
  vector<double> VMmult(const vector<int>& v, const aurostd::matrix<double>& A) {  // CO20200404 pflow::matrix()->aurostd::matrix()
    const uint M = A.size();
    const uint N = std::min((uint) A[0].size(), (uint) v.size());
    vector<double> u(M, 0.0);
    for (uint i = 0; i < N; i++) {
      double sum = 0;
      for (uint j = 0; j < M; j++) {
        sum = sum + v[j] * A[j][i];
      }
      u[i] = sum;
    }
    return u;
  }
} // namespace pflow

// ***************************************************************************
// Function VVcross
// ***************************************************************************
// This function returns the cross product of two 3 dimensional vectors. */
// Dane Morgan style
namespace pflow {
  vector<double> VVcross(const vector<double>& a, const vector<double>& b) {
    if (a.size() != 3 || b.size() != 3) {
      return a;
    }
    vector<double> u(3);
    u[0] = a[1] * b[2] - a[2] * b[1];
    u[1] = -(a[0] * b[2] - a[2] * b[0]);
    u[2] = a[0] * b[1] - a[1] * b[0];
    return u;
  }
} // namespace pflow

// ***************************************************************************
// Function VVdot
// ***************************************************************************
// This function returns the dot product of two vectors. */
// Dane Morgan style
namespace pflow {
  double VVdot(const vector<double>& a, const vector<double>& b) {
    double sum = 0;
    const int size = std::min(a.size(), b.size());
    for (int i = 0; i < size; i++) {
      sum = sum + a[i] * b[i];
    }
    return sum;
  }
} // namespace pflow

// ***************************************************************************
// Function GetNumAtoms
// ***************************************************************************
// This function returns the number of atoms in a structure
namespace pflow {
  int GetNumAtoms(const xstructure& a) {
    return a.atoms.size();
  }
} // namespace pflow

// ***************************************************************************
// Function SetSpline
// ***************************************************************************
// This function finds the second derivates so that
// spline interpolations can be evaluated quickly.
// I have shifted all indices to go from 0->(n-1).
namespace pflow {
  void SetSpline(const vector<double>& x, const vector<double>& y, const double& yp1, const double& ypn, vector<double>& y2) {
    int i;
    int k;
    double p;
    double qn;
    double sig;
    double un;
    const int n = x.size();

    vector<double> u(n - 1);
    y2 = vector<double>(n, 0.0);

    if (yp1 > 0.99e30) {
      y2[0] = u[0] = 0.0;
    } else {
      y2[0] = -0.5;
      u[0] = (3.0 / (x[2] - x[1])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }
    for (i = 1; i <= n - 2; i++) {
      sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
      p = sig * y2[i - 1] + 2.0;
      y2[i] = (sig - 1.0) / p;
      u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
      u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    if (ypn > 0.99e30) {
      qn = un = 0.0;
    } else {
      qn = 0.5;
      un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
    for (k = n - 2; k >= 0; k--) {
      y2[k] = y2[k] * y2[k + 1] + u[k];
    }
  }
} // namespace pflow

// ***************************************************************************
// Function GetSplineInt
// ***************************************************************************
// This function finds the value of a function y at a
// point x from a spline interpolation.  You must run
// SetSpline first.
// I have shifted all indices to go from 0->(n-1).
namespace pflow {
  void GetSplineInt(const vector<double>& xa, const vector<double>& ya, vector<double>& y2a, const double& x, double& y) {
    int klo;
    int khi;
    int k;
    double h;
    double b;
    double a;
    const int n = xa.size();

    klo = 0;
    khi = n - 1;
    while (khi - klo > 1) {
      k = (khi + klo) >> 1;
      if (xa[k] > x) {
        khi = k;
      } else {
        klo = k;
      }
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Bad xa input to routine splint", _INPUT_ERROR_);
    }
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
  }
} // namespace pflow

// ***************************************************************************
// Function PrintSpline
// ***************************************************************************
// This function prints out the original data and the
// spline interpolation with keywords for grep.
namespace pflow {
  void PrintSpline(const vector<double>& x, const vector<double>& y, const int& npts, ostream& outf) {
    outf.setf(std::ios::left, std::ios::adjustfield);
    outf.setf(std::ios::fixed, std::ios::floatfield);
    outf << "******************** Initial Data ********************" << endl;
    const int ndat = x.size();
    for (int id = 0; id < ndat; id++) {
      outf << x[id] << " " << y[id] << "   " << "INPUT" << endl;
    }
    outf << "******************** Cubic Spline Interpolated Points ********************" << endl;
    if (npts == 1) {
      outf << x[0] << " " << y[0] << "   " << "INTERPOLATED" << endl;
    } else {
      const double yp1 = 0.0;
      const double ypn = 0.0;
      vector<double> y2;
      SetSpline(x, y, yp1, ypn, y2);
      const double xrange = x[ndat - 1] - x[0];
      const double dx = xrange / (double) (npts - 1);
      for (int ip = 0; ip < npts; ip++) {
        const double xp = x[0] + (double) ip * dx;
        double yp;
        GetSplineInt(x, y, y2, xp, yp);
        outf << xp << " " << yp << "   " << "INTERPOLATED" << endl;
      }
    }
  }
} // namespace pflow

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
