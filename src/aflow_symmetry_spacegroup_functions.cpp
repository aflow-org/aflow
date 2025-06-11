// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Richard H. Taylor
// UPDATED by David Hicks
// d.hicks@duke.edu

#ifndef _AFLOW_SYMMETRY_SPACEGROUP_FUNCTIONS_CPP_
#define _AFLOW_SYMMETRY_SPACEGROUP_FUNCTIONS_CPP_
#include <string>
#include <vector>
#include <stdio.h>
#include <limits.h>  //Contains ULLONG_MAX value for precision check
#include <iostream>
#include <sstream>
#include "aflow.h"
#include "aflow_symmetry_spacegroup.h"
#include "aflow_compare_structure.h"
#include "aflow_pflow.h"
#include "aflow_contrib_shidong_cluster_expansion.h"

using namespace std;
using aurostd::FileExist;

//DX20190215 [OBSOLETE] #ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE] extern thread_local double _SYM_TOL_;
//DX20190215 [OBSOLETE] extern thread_local vector<xmatrix<double> > sym_mats;
//DX20190215 [OBSOLETE] extern thread_local vector<string> symbol;
//DX20190215 [OBSOLETE] #else
//DX20190215 [OBSOLETE] extern double _SYM_TOL_;
//DX20190215 [OBSOLETE] extern vector<xmatrix<double> > sym_mats;
//DX20190215 [OBSOLETE] extern vector<string> symbol;
//DX20190215 [OBSOLETE] #endif

//DX20190215 - added new class - START
// SymmetryInformationITC
// constructor
SymmetryInformationITC::SymmetryInformationITC(){
  //glides
  glideplanes.clear();
  glideplanes_hex.clear();
  glidetrans.clear();
  glidetrans_hex.clear();
  glidesymbols.clear();
  glidesymbols_hex.clear();
  //symmetry matrices
  index_cubic.clear();
  index_hex.clear();
  index_rhom.clear();
  index_tetr.clear();
  index_ortho.clear();
  index_mono_b.clear();
  index_mono_c.clear();
  index_tric.clear();
  sym_mats.clear();
  symbol.clear();
  dirparam.clear();
  sym_mats_direction.clear();
  sym_mats_hex.clear();
  symbol_hex.clear();
  dirparam_hex.clear();
  sym_mats_direction_hex.clear();
  //sym_ops
  sym_ops.clear();
  //generators
  generators.clear();
  sgindex.clear();
  gl_sgs.clear();
}

// free
void SymmetryInformationITC::free(){
}

// destructor
SymmetryInformationITC::~SymmetryInformationITC(){
  //glides
  glideplanes.clear();
  glideplanes_hex.clear();
  glidetrans.clear();
  glidetrans_hex.clear();
  glidesymbols.clear();
  glidesymbols_hex.clear();
  //symmetry matrices
  index_cubic.clear();
  index_hex.clear();
  index_rhom.clear();
  index_tetr.clear();
  index_ortho.clear();
  index_mono_b.clear();
  index_mono_c.clear();
  index_tric.clear();
  sym_mats.clear();
  symbol.clear();
  dirparam.clear();
  sym_mats_direction.clear();
  sym_mats_hex.clear();
  symbol_hex.clear();
  dirparam_hex.clear();
  sym_mats_direction_hex.clear();
  //sym_ops
  sym_ops.clear();
  //generators
  generators.clear();
  sgindex.clear();
  gl_sgs.clear();
  free();
}

// copy constructor
SymmetryInformationITC::SymmetryInformationITC(const SymmetryInformationITC& b){
  copy(b);
}

// copy constructor
void SymmetryInformationITC::copy(const SymmetryInformationITC& b){
  if(this != &b){
    glideplanes=b.glideplanes;
    glideplanes_hex=b.glideplanes_hex;
    glidetrans=b.glidetrans;
    glidetrans_hex=b.glidetrans_hex;
    glidesymbols=b.glidesymbols;
    glidesymbols_hex=b.glidesymbols_hex;
    //symmetry matrices
    index_cubic=b.index_cubic;
    index_hex=b.index_hex;
    index_rhom=b.index_rhom;
    index_tetr=b.index_tetr;
    index_ortho=b.index_ortho;
    index_mono_b=b.index_mono_b;
    index_mono_c=b.index_mono_c;
    index_tric=b.index_tric;
    sym_mats=b.sym_mats;
    symbol=b.symbol;
    dirparam=b.dirparam;
    sym_mats_direction=b.sym_mats_direction;
    sym_mats_hex=b.sym_mats_hex;
    symbol_hex=b.symbol_hex;
    dirparam_hex=b.dirparam_hex;
    sym_mats_direction_hex=b.sym_mats_direction_hex;
    //sym_ops
    sym_ops=b.sym_ops;
    //generators
    generators=b.generators;
    sgindex=b.sgindex;
    gl_sgs=b.gl_sgs;
  }
}

// ===== Assignment Operator (operator=) ===== //
const SymmetryInformationITC& SymmetryInformationITC::operator=(const SymmetryInformationITC& b){
  if(this!=&b){
    free();
    copy(b);
  }
  return *this;
}

// output
ostream& operator<<(ostream& oss, const SymmetryInformationITC& SymmetryInformationITC){

  oss << "SymmetryInformationITC: " << "Number of potential sym_mats: " << SymmetryInformationITC.sym_mats.size() << endl;
  oss << "SymmetryInformationITC: Cannot print all attributes; too much to print. Call attributes separately." << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "glideplanes: " << aurostd::joinWDelimiter(glideplanes,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "glideplanes_hex: " << aurostd::joinWDelimiter(glideplanes_hex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "glidetrans: " << aurostd::joinWDelimiter(glidetrans,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "glidetrans_hex: " << aurostd::joinWDelimiter(glidetrans_hex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "glidesymbols: " << aurostd::joinWDelimiter(glidesymbols,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "glidesymbols_hex: " << aurostd::joinWDelimiter(glidesymbols_hex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "index_cubic: " << aurostd::joinWDelimiter(index_cubic,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "index_hex: " << aurostd::joinWDelimiter(index_hex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "index_rhom: " << aurostd::joinWDelimiter(index_rhom,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "index_tetr: " << aurostd::joinWDelimiter(index_tetr,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "index_ortho: " << aurostd::joinWDelimiter(index_ortho,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "index_mono_b: " << aurostd::joinWDelimiter(index_mono_b,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "index_mono_c: " << aurostd::joinWDelimiter(index_mono_c,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "index_tric: " << aurostd::joinWDelimiter(index_tric,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "sym_mats: " << aurostd::joinWDelimiter(sym_mats,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "symbol: " << aurostd::joinWDelimiter(symbol,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "dirparam: " << aurostd::joinWDelimiter(dirparam,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "sym_mats_direction: " << aurostd::joinWDelimiter(sym_mats_direction,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "sym_mats_hex: " << aurostd::joinWDelimiter(sym_mats_hex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "symbol_hex: " << aurostd::joinWDelimiter(symbol_hex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "dirparam_hex: " << aurostd::joinWDelimiter(dirparam_hex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "sym_mats_direction_hex: " << aurostd::joinWDelimiter(sym_mats_direction_hex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "sym_ops: " << aurostd::joinWDelimiter(sym_ops,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "generators: " << aurostd::joinWDelimiter(generators,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "sgindex: " << aurostd::joinWDelimiter(sgindex,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]oss << "gl_sgs: " << aurostd::joinWDelimiter(gl_sgs,",") << endl;
  //DX20190215 [TOO MUCH TO PRINT]
  return oss;
}
//DX20190215 - added new class - END


// ******************************************************************************
// getPearsonSymbol
// ******************************************************************************
// Determine the Pearson symbol of the newly found conventional cell
namespace SYM {
  string getPearsonSymbol(char& centering, char& lattice_char, deque<_atom> atoms) {
    ostringstream psymb;
    char pearson_lattice_char = lattice_char;
    // === Change lattice_char if monoclinic unique axis b was used to 'm' === //
    if(lattice_char == 'b') {
      pearson_lattice_char = 'm';
    }
    // === Pearson symbol centering uses 'S' for "side face centered" (i.e. S=A,B,or C) === //
    char pearson_centering = centering;
    if(centering == 'A' || centering == 'B' || centering == 'C') {
      pearson_centering = 'S';
    }
    psymb << pearson_lattice_char << pearson_centering << atoms.size();
    return psymb.str();
  }
} //namespace SYM

//DX20190418 - START
// ******************************************************************************
// spacegroup2latticeAndCentering
// ******************************************************************************
// Categorize space group into the corresponding lattice and centering (1 of 14)
namespace SYM {
  string spacegroup2latticeAndCentering(uint sg_number){

    // ---------------------------------------------------------------------------
    // triclinic
    if(sg_number >= 1 && sg_number <= 2){
      return "aP";
    }

    // ---------------------------------------------------------------------------
    // monoclinic
    if(sg_number >= 3 && sg_number <= 15){
      // simple monoclinic
      if(sg_number == 3 || sg_number == 4 || sg_number == 6 || sg_number == 7 ||
          sg_number == 10 || sg_number == 11 || sg_number == 13 || sg_number == 14){
        return "mP";
      }
      // base-centered monoclinic
      else if(sg_number == 5 || sg_number == 8 || sg_number == 9 || sg_number == 12 || sg_number == 15){
        return "mC";
      }
    }

    // ---------------------------------------------------------------------------
    // orthrohombic
    if(sg_number >= 16 && sg_number <= 74){
      // simple orthorhombic
      if ((sg_number >= 16 && sg_number <= 19) || (sg_number >= 25 && sg_number <= 34) || (sg_number >= 47 && sg_number <= 62)){
        return "oP";
      }
      // base-centered orthorhombic
      else if((sg_number >= 20 && sg_number <= 21) || (sg_number >= 35 && sg_number <= 41) || (sg_number >= 63 && sg_number <= 68)){
        return "oC";
      }
      // face-centered orthorhombic
      else if((sg_number == 22) || (sg_number >= 42 && sg_number <= 43) || (sg_number >= 69 && sg_number <= 70)){
        return "oF";
      }
      // body-centered orthorhombic
      else if((sg_number >= 23 && sg_number <= 24) || (sg_number >= 44 && sg_number <= 46) || (sg_number >= 71 && sg_number <= 74)){
        return "oI";
      }
    }

    // ---------------------------------------------------------------------------
    // tetragonal
    if(sg_number >= 75 && sg_number <= 142){
      // simple tetragonal
      if((sg_number >= 75 && sg_number <= 78) || (sg_number == 81) || (sg_number >= 83 && sg_number <= 86) ||
          (sg_number >= 89 && sg_number <= 96) || (sg_number >= 99 && sg_number <= 106) || (sg_number >= 111 && sg_number <= 118) ||
          (sg_number >= 123 && sg_number <= 138)){
        return "tP";
      }
      // body-centered tetragonal
      else if((sg_number >= 79 && sg_number <= 80) || (sg_number == 82) || (sg_number >= 87 && sg_number <= 88) ||
          (sg_number >= 97 && sg_number <= 98) || (sg_number >= 107 && sg_number <= 110) || (sg_number >= 119 && sg_number <= 122) ||
          (sg_number >= 139 && sg_number <= 142)){
        return "tI";
      }
    }

    // ---------------------------------------------------------------------------
    // trigonal
    if(sg_number >= 143 && sg_number <= 167){
      // trigonal
      if((sg_number >= 143 && sg_number <= 145) || (sg_number == 147) || (sg_number >= 149 && sg_number <= 154) ||
          (sg_number >= 156 && sg_number <= 159) || (sg_number >= 162 && sg_number <= 165)){
        return "hP";
      }
      // rhombohedral
      else if((sg_number == 146) || (sg_number == 148) || (sg_number == 155) ||
          (sg_number >= 160 && sg_number <= 161) || (sg_number >= 166 && sg_number <= 167)){
        return "hR";
      }
    }

    // ---------------------------------------------------------------------------
    // hexagonal
    if(sg_number >= 168 && sg_number <= 194){
      return "hP";
    }

    // ---------------------------------------------------------------------------
    // cubic
    if(sg_number >= 195 && sg_number <= 230){
      // simple cubic
      if(sg_number == 195 || sg_number == 198 || sg_number == 200 || sg_number == 201 || sg_number == 205 || sg_number == 207 ||
          sg_number == 208 || sg_number == 212 || sg_number == 213 || sg_number == 215 || sg_number == 218 || sg_number == 221 ||
          sg_number == 222 || sg_number == 223 || sg_number == 224){
        return "cP";
      }
      // face-centered cubic
      else if(sg_number == 196 || sg_number == 202 || sg_number == 203 || sg_number == 209 || sg_number == 210 || sg_number == 216 ||
          sg_number == 219 || sg_number == 225 || sg_number == 226 || sg_number == 227 || sg_number == 228){
        return "cF";
      }
      // body-centered cubic
      else if(sg_number == 197 || sg_number == 199 || sg_number == 204 || sg_number == 206 || sg_number == 211 || sg_number == 214 ||
          sg_number == 217 || sg_number == 220 || sg_number == 229 || sg_number == 230){
        return "cI";
      }
    }

    // ---------------------------------------------------------------------------
    // error
    stringstream message;
    message << "Could not classify space group into a Bravais lattice (most likely ranges are incorrect).";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);

  }
}
//DX20190418 - END

// ******************************************************************************
// isEnantimorphicPair
// ******************************************************************************
// Return enantiomorph of space group.
// If none, exists, return original space group
namespace SYM {
  uint getEnantiomorphSpaceGroupNumber(uint space_group_number){
    if(space_group_number==76){return 78;}
    if(space_group_number==78){return 76;}
    if(space_group_number==91){return 95;}
    if(space_group_number==95){return 91;}
    if(space_group_number==92){return 96;}
    if(space_group_number==96){return 92;}
    if(space_group_number==144){return 145;}
    if(space_group_number==145){return 144;}
    if(space_group_number==151){return 153;}
    if(space_group_number==153){return 151;}
    if(space_group_number==152){return 154;}
    if(space_group_number==154){return 152;}
    if(space_group_number==169){return 170;}
    if(space_group_number==170){return 169;}
    if(space_group_number==178){return 179;}
    if(space_group_number==179){return 178;}
    if(space_group_number==180){return 181;}
    if(space_group_number==181){return 180;}
    if(space_group_number==212){return 213;}
    if(space_group_number==213){return 212;}
    return space_group_number;
  }
}

// ******************************************************************************
// getAtomGCD
// ******************************************************************************
// Determine the greatest common denominator (GCD) between the number of atoms in
// the primitive atomic basis. This is used to check the consistency between the
// ratio of atoms when we find a conventional cell
namespace SYM {
  bool getAtomGCD(deque<_atom>& atomic_basis, deque<deque<_atom> >& split_atom_types, int& gcd_num) {
    split_atom_types = SYM::break_up_by_type(atomic_basis);
    if(split_atom_types.size()==0){throw aurostd::xerror(__AFLOW_FILE__,"SYM::getAtomGCD():","split_atom_types.size()==0",_INPUT_ERROR_);}
    gcd_num=1;  //CO+DX20191201
    if(split_atom_types.size() > 1) { // If only one type of atom, set to 1
      gcd_num=(int)split_atom_types[0].size();  //CO+DX20191201
      for (uint p = 0; p < split_atom_types.size() - 1; p++) {
        aurostd::GCD(gcd_num, (int)split_atom_types[p + 1].size(),gcd_num);  //CO20191201
      }
    }
    return true;
  }
}

// ******************************************************************************
// latticeVectorsSame
// ******************************************************************************
// Check if the lattice vectors are the same
namespace SYM{
  bool latticeVectorsSame(xvector<double>& a, xvector<double>& b, xvector<double>& c, double& tol) {
    return ((aurostd::abs(a(1) - b(1)) <= tol && aurostd::abs(a(2) - b(2)) <= tol && aurostd::abs(a(3) - b(3)) <= tol) ||
        (aurostd::abs(a(1) - c(1)) <= tol && aurostd::abs(a(2) - c(2)) <= tol && aurostd::abs(a(3) - c(3)) <= tol) ||
        (aurostd::abs(b(1) - c(1)) <= tol && aurostd::abs(b(2) - c(2)) <= tol && aurostd::abs(b(3) - c(3)) <= tol));
  }
} //namespace SYM

// ******************************************************************************
// orientVectorsPositiveQuadrant
// ******************************************************************************
// Orient vectors in the positive quadrant
namespace SYM {
  bool orientVectorsPositiveQuadrant(vector<xvector<double> > lattice_vectors, double& tol) {
    for (uint i = 0; i < lattice_vectors.size(); i++) {
      orientVectorsPositiveQuadrant(lattice_vectors[i], tol);
    }
    return true;
  }
} //namespace SYM

namespace SYM {
  bool orientVectorsPositiveQuadrant(xvector<double>& vec, double& tol) {
    if(vec(1) < -tol || vec(2) < -tol || vec(3) < -tol) {
      vec = -1.0 * vec;
    }
    return true;
  }
} //namespace SYM

// ******************************************************************************
// alignLatticeWithXYZ
// ******************************************************************************
// Align lattice with XYZ coordinates
namespace SYM {
  bool alignLatticeWithXYZ(xvector<double>& a, xvector<double>& b, xvector<double>& c, double& tol) {
    if(b(1) > tol && b(2) < tol && b(3) < tol) {
      xvector<double> tmp = a;
      a = b;
      b = tmp;
    } else if(c(1) > tol && c(2) < tol && c(3) < tol) {
      xvector<double> tmp = a;
      a = c;
      c = tmp;
    }
    if(a(1) < tol && a(2) > tol && a(3) < tol) {
      xvector<double> tmp = b;
      b = a;
      a = tmp;
    } else if(c(1) < tol && c(2) > tol && c(3) < tol) {
      xvector<double> tmp = b;
      b = c;
      c = tmp;
    }
    return true;
  }
} //namespace SYM

// ******************************************************************************
// anyNullVectors
// ******************************************************************************
// Check for any null vectors
namespace SYM {
  bool anyNullVectors(vector<xvector<double> >& vecs, double& tol) {
    for (uint v = 0; v < vecs.size(); v++) {
      if(nullVector(vecs[v], tol)) {
        return true;
      }
    }
    return false;
  }
} //namespace SYM

namespace SYM {
  bool nullVector(xvector<double>& vec, double& tol) {
    return (aurostd::abs(vec[1]) < tol && aurostd::abs(vec[2]) < tol && aurostd::abs(vec[3]) < tol);
  }
} //namespace SYM

// ******************************************************************************
// updateAtomPositions
// ******************************************************************************
// Update atom positions after screw operation
namespace SYM {
  void updateAtomPositions(deque<_atom>& atoms, Screw& S, xmatrix<double>& lattice) { //DX20190905 - return to void 
    for (uint a = 0; a < atoms.size(); a++) {
      atoms[a].cpos = S * atoms[a].cpos;
      //DX20190905 [OBSOLETE-no more mod_one_xvec] atoms[a].fpos = SYM::mod_one_xvec(C2F(lattice, atoms[a].cpos));
      atoms[a].fpos = C2F(lattice, atoms[a].cpos); //DX20190905
      BringInCellInPlace(atoms[a].fpos); //DX20190905
    }
    //DX20190905 [OBSOLETE] return atoms;
  }
} //namespace SYM

// ******************************************************************************
// Dir
// ******************************************************************************
namespace SYM {
  xvector<double> dir(double a, double b, double c) {
    xvector<double> out;
    out(1) = a;
    out(2) = b;
    out(3) = c;
    return out;
  }
}

//DX20190128 START
// ******************************************************************************
// ExtractLatticeParametersFromWyccar
// ******************************************************************************
namespace SYM {
  vector<double> ExtractLatticeParametersFromWyccar(const vector<string>& wyccar_ITC){ //DX20191030 - added const
    vector<double> lattice_parameters;
    if(wyccar_ITC.size()>3){
      vector<string> tokens;
      aurostd::string2tokens(wyccar_ITC[2],tokens," ");
      for(uint i=0;i<tokens.size();i++){
        if(i<6){ lattice_parameters.push_back(aurostd::string2utype<double>(tokens[i])); }
      }
    }
    else {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid WYCCAR",_FILE_WRONG_FORMAT_);
    }
    return lattice_parameters;
  }
}
//DX20190128 END

// ******************************************************************************
// ExtractWyckoffAttributeString 
// ******************************************************************************
namespace SYM {
  string ExtractWyckoffAttributeString(const vector<string>& wyccar_ITC, uint attribute_index){ //DX20191030 - added const
    vector<string> all_wyckoff_sets;
    vector<string> wyckoff_set;
    string atom_name = ""; string prev_atom_name = ""; string wyckoff_letter = "";
    for(uint i=0;i<wyccar_ITC.size();i++){
      if(i>4 && i!=wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
        vector<string> tokens;
        aurostd::string2tokens(wyccar_ITC[i],tokens," ");
        string atom_name = tokens[3];
        wyckoff_letter = tokens[attribute_index];
        if(atom_name != prev_atom_name){
          if(wyckoff_set.size()>0){
            all_wyckoff_sets.push_back(aurostd::joinWDelimiter(wyckoff_set,","));
            wyckoff_set.clear(); 
          }
          prev_atom_name = atom_name;  
        }
        wyckoff_set.push_back(wyckoff_letter);
      }
    }
    all_wyckoff_sets.push_back(aurostd::joinWDelimiter(wyckoff_set,","));
    wyckoff_set.clear(); 
    return aurostd::joinWDelimiter(all_wyckoff_sets,";");    
  }
}

// ******************************************************************************
// ExtractWyckoffLettersString 
// ******************************************************************************
namespace SYM {
  string ExtractWyckoffLettersString(const vector<string>& wyccar_ITC){ //DX20191030 - added const
    return ExtractWyckoffAttributeString(wyccar_ITC, 5); //Wyckoff letter sits at 5th column of wyccar
  }
}

// ******************************************************************************
// ExtractWyckoffMultiplicitiesString 
// ******************************************************************************
namespace SYM {
  string ExtractWyckoffMultiplicitiesString(const vector<string>& wyccar_ITC){ //DX20191030 - added const
    return ExtractWyckoffAttributeString(wyccar_ITC, 4); //Wyckoff multiplicity sits at 4th column of wyccar
  }
}

// ******************************************************************************
// ExtractWyckoffSiteSymmetriesString 
// ******************************************************************************
namespace SYM {
  string ExtractWyckoffSiteSymmetriesString(const vector<string>& wyccar_ITC){ //DX20191030 - added const
    return ExtractWyckoffAttributeString(wyccar_ITC, 6); //Wyckoff site symmetry sits at 6th column of wyccar
  }
}

//DX20190128 START
// ******************************************************************************
// getWyckoffLettersWithSameMultiplcityAndSiteSymmetry 
// ******************************************************************************
namespace SYM {
  vector<vector<vector<string> > > getWyckoffLettersWithSameMultiplcityAndSiteSymmetry(uint& space_group_number, 
      vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, uint& cell_choice){ 
    // Get the set of Wyckoff letters with the same multiplicity and site symmetry
    bool LDEBUG=(false || XHOST.DEBUG);

    stringstream axis_cell;
    axis_cell.str(std::string());
    axis_cell << cell_choice; 
    if(LDEBUG) { cerr << "SYM::getWyckoffLettersWithSameMultiplcityAndSiteSymmetry(): Lattice type (unique axis if specified) and cell choice  : " << axis_cell.str() << endl; }
    SymmetryInformationITC ITC_sym_info; //DX20190215
    ITC_sym_info.initsgs(axis_cell.str()); //DX20190215

    //DX20190215 [OBSOLETE] using SYM::gl_sgs;
    string spacegroupstring = ITC_sym_info.gl_sgs[space_group_number - 1];
    vector<vector<vector<string> > > grouped_possible_Wyckoff_letters;
    for(uint i=0;i<grouped_Wyckoff_positions.size();i++){
      vector<vector<string> > possible_Wyckoff_letters, grouped_possible_Wyckoff_site_symmetries;
      for(uint j=0;j<grouped_Wyckoff_positions[i].site_symmetries.size();j++){
        vector<string> site_syms, letters, positions;
        SYM::get_certain_wyckoff_pos(spacegroupstring, grouped_Wyckoff_positions[i].multiplicities[j], grouped_Wyckoff_positions[i].site_symmetries[j], site_syms, letters, positions);
        possible_Wyckoff_letters.push_back(letters);
      }
      grouped_possible_Wyckoff_letters.push_back(possible_Wyckoff_letters);
    }

    // print
    if(LDEBUG) {
      for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){
        cerr << "species " << i << ": " << endl;;
        for(uint j=0;j<grouped_possible_Wyckoff_letters[i].size();j++){
          cerr << "  possible Wyckoff letters for group " << i << ": " << aurostd::joinWDelimiter(grouped_possible_Wyckoff_letters[i][j],",") << endl;;
        }
      }
    }

    return grouped_possible_Wyckoff_letters;
  }
}
//DX20190128 END

//DX20190219 START
// ******************************************************************************
// splitSiteSymmetry 
// ******************************************************************************
namespace SYM {
  vector<string> splitSiteSymmetry(const string& site_symmetry){ //DX20190730 - added const

    // Split the site symmetry symbol into its primary, secondary, and teriary
    // directions (see ITC-A pgs. 28-29 for site symmetry info).


    vector<string> split_site_symmetry; 
    stringstream ss_primary, ss_secondary, ss_tertiary;

    // if only primary direction has symmetry
    if(site_symmetry.size()==1){
      split_site_symmetry.push_back(site_symmetry);
      return split_site_symmetry;
    }

    // if only primary direction (-1) or primary and secondary only (happens for rhl)
    if(site_symmetry.size()==2){
      // if inversion (-1)
      if(site_symmetry[0]=='-'){
        split_site_symmetry.push_back(site_symmetry);
        return split_site_symmetry;
      }
      else {
        ss_primary << site_symmetry[0]; ss_secondary << site_symmetry[1];
        split_site_symmetry.push_back(ss_primary.str());
        split_site_symmetry.push_back(ss_secondary.str());
        ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
        return split_site_symmetry;
      }
    }

    // if only primary direction with n-fold with perpendicular mirror
    if(site_symmetry.size()==3 && aurostd::substring2bool(site_symmetry,"/m")){
      split_site_symmetry.push_back(site_symmetry);
      return split_site_symmetry;
    }

    // if only primary and secondary direction with inversion 
    if(site_symmetry.size()==3 && aurostd::substring2bool(site_symmetry,"-")){
      if(site_symmetry[0]=='-'){
        ss_primary << site_symmetry[0] << site_symmetry[1]; ss_secondary << site_symmetry[2];
        split_site_symmetry.push_back(ss_primary.str());
        split_site_symmetry.push_back(ss_secondary.str());
        ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
        return split_site_symmetry;
      }
      else {
        ss_primary << site_symmetry[0]; ss_secondary << site_symmetry[1] << site_symmetry[2];
        split_site_symmetry.push_back(ss_primary.str());
        split_site_symmetry.push_back(ss_secondary.str());
        ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
        return split_site_symmetry;
      }
    }

    // check non-symmetry directions
    vector<uint> dot_index;
    vector<string> dot_direction; // primary, secondary, tertiary
    if(aurostd::substring2bool(site_symmetry,".")){
      for(uint i=0;i<site_symmetry.size();i++){
        if(site_symmetry[i] == '.'){
          dot_index.push_back(i);
          if(i==0){ dot_direction.push_back("primary"); }
          else if(i<site_symmetry.size()-1){ 
            dot_direction.push_back("secondary");
          }
          else if(i==site_symmetry.size()-1){
            dot_direction.push_back("tertiary");
          }
        }
      }

      // if two non-symmetry directions
      if(dot_index.size()==2){
        if(dot_direction[0]=="primary"){ 
          ss_primary << site_symmetry[dot_index[0]];
          split_site_symmetry.push_back(ss_primary.str());
          ss_primary.str("");
          if(dot_direction[1]=="secondary"){
            // form: ..<sym>
            ss_secondary << site_symmetry[dot_index[1]];
            ss_tertiary << site_symmetry.substr(dot_index[1]+1,site_symmetry.size());
            split_site_symmetry.push_back(ss_secondary.str());
            split_site_symmetry.push_back(ss_tertiary.str());
            ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
          }
          else if(dot_direction[1]=="tertiary"){
            // form: .<sym>.
            ss_secondary << site_symmetry.substr(dot_index[0]+1,dot_index[1]-1);
            ss_tertiary << site_symmetry[dot_index[1]];
            split_site_symmetry.push_back(ss_secondary.str());
            split_site_symmetry.push_back(ss_tertiary.str());
            ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
          }
        }
        else {
          // form: <sym>..
          ss_primary << site_symmetry.substr(0,dot_index[0]);
          ss_secondary << site_symmetry[dot_index[0]];
          ss_tertiary << site_symmetry[dot_index[1]];
          split_site_symmetry.push_back(ss_primary.str());
          split_site_symmetry.push_back(ss_secondary.str());
          split_site_symmetry.push_back(ss_tertiary.str());
          ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
        }
      }

      // if one non-symmetry direction
      if(dot_index.size()==1){
        if(dot_direction[0]=="primary"){ 
          ss_primary << site_symmetry[dot_index[0]];
          split_site_symmetry.push_back(ss_primary.str());
          ss_primary.str("");
          string tmp = site_symmetry.substr(dot_index[0]+1,site_symmetry.size());
          if(tmp.size()==2){
            ss_secondary << tmp[0];
            ss_tertiary << tmp[1];
            split_site_symmetry.push_back(ss_secondary.str());
            split_site_symmetry.push_back(ss_tertiary.str());
            ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
          }
          else if(tmp.size()==3){
            if(isdigit(tmp[0])){ //ex .2mm
              ss_secondary << tmp[0] << tmp[1];
              ss_tertiary << tmp[2];
              split_site_symmetry.push_back(ss_secondary.str());
              split_site_symmetry.push_back(ss_tertiary.str());
              ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
            }
            else if(!isdigit(tmp[0])){ //ex .m2m
              ss_secondary << tmp[0];
              ss_tertiary << tmp[1] << tmp[2];
              split_site_symmetry.push_back(ss_secondary.str());
              split_site_symmetry.push_back(ss_tertiary.str());
              ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
            }
            else { // not sure, but by default, higher sym goes to secondary direction
              ss_secondary << tmp[0] << tmp[1];
              ss_tertiary << tmp[2];
              split_site_symmetry.push_back(ss_secondary.str());
              split_site_symmetry.push_back(ss_tertiary.str());
              ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
            }
          }
          else if(tmp.size()==4){
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"More site symmetry characters than anticipated.",_INPUT_ILLEGAL_);
          }
        }
        if(dot_direction[0]=="secondary"){ 
          ss_primary << site_symmetry.substr(0,dot_index[0]);
          ss_secondary << site_symmetry[dot_index[0]];
          ss_tertiary << site_symmetry.substr(dot_index[0]+1,site_symmetry.size());
          split_site_symmetry.push_back(ss_primary.str());
          split_site_symmetry.push_back(ss_secondary.str());
          split_site_symmetry.push_back(ss_tertiary.str());
          ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
        }
        if(dot_direction[0]=="tertiary"){ 
          string tmp = site_symmetry.substr(0,dot_index[0]);
          if(tmp.size()==2){
            ss_primary << tmp[0];
            ss_secondary << tmp[1];
            split_site_symmetry.push_back(ss_primary.str());
            split_site_symmetry.push_back(ss_secondary.str());
            ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
          }
          else if(tmp.size()==3){
            if(isdigit(tmp[0])){ //ex 2mm.
              ss_primary << tmp[0];
              ss_secondary << tmp[1] << tmp[2];
              split_site_symmetry.push_back(ss_primary.str());
              split_site_symmetry.push_back(ss_secondary.str());
              ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
            }
            else if(!isdigit(tmp[0])){ //ex m2m.
              ss_primary << tmp[0];
              ss_secondary << tmp[1] << tmp[2];
              split_site_symmetry.push_back(ss_primary.str());
              split_site_symmetry.push_back(ss_secondary.str());
              ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
            }
            else { 
              ss_primary << tmp[0];
              ss_secondary << tmp[1] << tmp[2];
              split_site_symmetry.push_back(ss_primary.str());
              split_site_symmetry.push_back(ss_secondary.str());
              ss_primary.str(""); ss_secondary.str(""); ss_tertiary.str("");
            }
          }
          else if(tmp.size()==4){
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"More site symmetry characters than anticipated.",_INPUT_ILLEGAL_);
          }
          ss_tertiary << site_symmetry[dot_index[0]];
          split_site_symmetry.push_back(ss_tertiary.str());
          ss_tertiary.str("");
        }
      }
    }
    else {
      bool is_inversion = false;
      bool perp_mirror = false;
      stringstream ss_tmp;
      for(uint i=0;i<site_symmetry.size();i++){
        if(site_symmetry[i]=='-'){
          ss_tmp << site_symmetry[i];
          is_inversion=true;
        }
        else if(is_inversion){
          ss_tmp << site_symmetry[i];
          is_inversion=false;
          split_site_symmetry.push_back(ss_tmp.str());
          ss_tmp.str("");
        }
        else if(site_symmetry[i]=='/'){
          ss_tmp << site_symmetry[i];
          perp_mirror=true;
        }
        else if(perp_mirror){
          ss_tmp << site_symmetry[i];
          perp_mirror=false;
          split_site_symmetry.push_back(ss_tmp.str());
          ss_tmp.str("");
        }
        else {
          ss_tmp << site_symmetry[i];
          if(i<site_symmetry.size()-1){
            if(site_symmetry[i+1]=='/'){
              continue;
            }
          }
          split_site_symmetry.push_back(ss_tmp.str());
          ss_tmp.str("");
        }
      }
    }

    if(split_site_symmetry.size()!=3){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"More site symmetry characters than anticipated.",_INPUT_ILLEGAL_);
    }

    return split_site_symmetry;
  }
}


//// ******************************************************************************
//// Atomic Envirionment
//// ******************************************************************************
//void AtomicEnvironment(istream & cin, vector<string> av ){
//
//  if(av.size() != 3) {throw aurostd::xerror(__AFLOW_FILE__,XPID+"AtomicEnvironment():","Please enter a radius.",_INPUT_NUMBER_);}
//  double radius = atof(av[2].c_str());
//  //double tol = 1e-6;
//  CrystalStructure C(cin);//Gives primitive lattice
//  xmatrix<double> L = C.return_prim_lattice();
//  xmatrix<double> Linv = aurostd::inverse(L);
//
//  vector<xvector<double> > expanded_lattice;
//  vector<_atom> expanded_basis;
//  //Scale the expansion to include all the points up to the radius given by user
//  int b=0;
//  //double scale;
////DX  double lp = C.return_lattice_parameter();
//  double lp = C.return_scale();
//  double min;
//  //If volume:
////DX HERE  if(lp < 0){
////DX HERE   min = modulus(extract_row(L,1));
////DX HERE   if(modulus(extract_row(L,2)) < min){
////DX HERE     min = modulus(extract_row(L,2));
////DX HERE   }
////DX HERE   if(modulus(extract_row(L,3)) < min) {
////DX HERE     min = modulus(extract_row(L,3));
////DX HERE   }
////DX HERE    b = (int) (radius/min+1.0);
////DX HERE  }
//  //If lattice parameter:
////DX HERE  if(lp > 0){
//    min = modulus(lp*extract_row(L,1));
//    if(modulus(lp*extract_row(L,2)) < min){
//      min = modulus(lp*extract_row(L,2));
//    }
//    if(modulus(lp*extract_row(L,3)) < min) {
//      min = modulus(lp*extract_row(L,3));
//    }
//    b = (int) (radius/min+1.0);
////DX HERE  }
//  //EXPAND BASIS:
//  expanded_lattice = C.expand_lattice(b,b,b,L);
//  expanded_basis = C.add_basis(expanded_lattice);
//
////FIND DIFERENCES BETWEEN EXPANSION AND ATOMS IN PRIMITIVE
//  vector<_atom> atomic_basis_ = C.return_atomic_basis_();
//  //CALCULATE DISTANCES
//  vector<vector<stringdouble> > all_dist_vec;
//  vector<stringdouble> one_dist_vec;
//  stringdouble one_dist;
//  for(uint i=0;i<atomic_basis_.size();i++){
//    for(uint j=0;j<expanded_basis.size();j++){
//      one_dist.distance = distance_between_points(expanded_basis[j].cpos, atomic_basis_[i].fpos*L);
//      one_dist.type = expanded_basis[j].type;
//      one_dist.coord = expanded_basis[j].coord;
//      one_dist_vec.push_back(one_dist);
//    }
//    all_dist_vec.push_back(one_dist_vec);
//    one_dist_vec.clear();
//  }
//
//
//////////////////////////////////////////////////////////////////////
//  //bubble sort
//
//  for(uint i=0;i<all_dist_vec.size();i++){
//    for(uint j=0;j<all_dist_vec[i].size();j++){
//      for(uint k=0;k<all_dist_vec[i].size();k++){
//	//If the kth atom is closer than the ith atom switch them:
//      	if(all_dist_vec[i][j].distance - all_dist_vec[i][k].distance < 0){
//	  stringdouble tmp;
//	  tmp = all_dist_vec[i][k];
//	  all_dist_vec[i][k] = all_dist_vec[i][j];
//	  all_dist_vec[i][j] = tmp;
//	}
//      }
//    }
//  }
//  for(uint i=0;i<all_dist_vec.size();i++){
//    cerr << "*******************" << endl;
//    cerr << "ATOM " << i << endl;
//    cerr << "*******************" << endl;
//    for(uint j=0;j<all_dist_vec[i].size();j++){
//      if(all_dist_vec[i][j].distance < radius){
//	cerr << all_dist_vec[i][j].distance << endl;
//      }
//      //      cout <<"relative distance (wrt nn): " <<  all_dist_vec[i][j].distance/all_dist_vec[i][1].distance << " " << all_dist_vec[i][j].type << endl;
//    }
//  }
//}

// ******************************************************************************
// Orthogonality Defect
// ******************************************************************************
// Calculate the orthogonality defect
namespace SYM {
  double orthogonality_defect(xmatrix<double> lattice) {
    double a = sqrt(lattice(1, 1) * lattice(1, 1) + lattice(1, 2) * lattice(1, 2) + lattice(1, 3) * lattice(1, 3));
    double b = sqrt(lattice(2, 1) * lattice(2, 1) + lattice(2, 2) * lattice(2, 2) + lattice(2, 3) * lattice(2, 3));
    double c = sqrt(lattice(3, 1) * lattice(3, 1) + lattice(3, 2) * lattice(3, 2) + lattice(3, 3) * lattice(3, 3));
    double det = aurostd::abs(aurostd::determinant(lattice));
    return a * b * c / det;
  }
} //namespace SYM

// ******************************************************************************
// Length Lattice Vectors
// ******************************************************************************
// Sum of lattice vectors lengths
namespace SYM {
  double length_lattice_vectors(xmatrix<double> lattice) {
    double a = sqrt(lattice(1, 1) * lattice(1, 1) + lattice(1, 2) * lattice(1, 2) + lattice(1, 3) * lattice(1, 3));
    double b = sqrt(lattice(2, 1) * lattice(2, 1) + lattice(2, 2) * lattice(2, 2) + lattice(2, 3) * lattice(2, 3));
    double c = sqrt(lattice(3, 1) * lattice(3, 1) + lattice(3, 2) * lattice(3, 2) + lattice(3, 3) * lattice(3, 3));
    return a + b + c;
  }
} //namespace SYM

// ******************************************************************************
// Smallest, Greater than Minimum
// ******************************************************************************
// Smallest number in vec that is greater than min
namespace SYM {
  double smallest_gt_min(double min, vector<double> vec) {
    double eps = 1e-6;
    double out = aurostd::max(vec); //DX20210422 - use sym::infnorm -> aurostd::max
    for (uint i = 0; i < vec.size(); i++) {
      if((vec[i] - min) > eps && vec[i] <= out) {
        out = vec[i];
      }
    }
    return out;
  }
} //namespace SYM

// ******************************************************************************
// Smallest, Greater than Minimum Indicies
// ******************************************************************************
// Smallest number in vec that is greater than min
namespace SYM {
  int smallest_gt_min_index(double min, int not_index1, int not_index2, vector<double> vec) {
    double eps = 1e-6;
    double out = aurostd::max(vec); //DX20210422 - use sym::infnorm -> aurostd::max
    int index = 0; //DX20180514 - added initialization 
    for (int i = 0; (uint)i < vec.size(); i++) {
      if((vec[i] - min) > eps && i != not_index1 && i != not_index2 && vec[i] <= out) {
        out = vec[i];
        index = i;
      } else if(aurostd::abs(vec[i] - min) <= eps && i != not_index1 && i != not_index2) {
        out = vec[i];
        index = i;
      }
    }
    return index;
  }
} //namespace SYM

// ******************************************************************************
// Cleanup String
// ******************************************************************************
namespace SYM {
  void cleanupstring(string& str) {
    if(str.size() > 0) {
      uint i = str.size() - 1;
      while (str[i] == ' ') {
        i--;
      }
      str.erase(i + 1, str.size() - 1);
      i = 0;
      while (str[i] == ' ') {
        i++;
      }
      str.erase(0, i);
    }
  }
} //namespace SYM

// ******************************************************************************
// Split String By Spaces
// ******************************************************************************
namespace SYM {
  vector<string> splitstringbyspaces(string str) {
    str.append(" ");
    vector<string> out;
    ostringstream oss;
    for (uint i = 0; i < str.size(); i += 0) {
      if(str[i] != ' ') {
        oss << str[i];
      }
      if(i > 0 && str[i] == ' ') {
        out.push_back(oss.str());
        oss.str("");
        while (str[i] == ' ') {
          i++;
        }
      } else {
        i++;
      }
    }
    return out;
  }
} //namespace SYM

// ******************************************************************************
// Split String
// ******************************************************************************
//Split a string with elements separated by delimeter 'c.' Antiquates vector<string> splitstringbyspaces.
namespace SYM {
  vector<string> splitstring(string str, char c) {
    str += c;
    vector<string> out;
    ostringstream oss;
    for (uint i = 0; i < str.size(); i += 0) {
      if(str[i] != c) {
        oss << str[i];
      }
      if(i > 0 && str[i] == c) {
        out.push_back(oss.str());
        oss.str("");
        while (str[i] == c) {
          i++;
        }
      } else {
        i++;
      }
    }
    return out;
  }
} //namespace SYM

// ******************************************************************************
// ReducedRowEchelonForm
// ******************************************************************************
namespace SYM {
  void ReducedRowEchelonForm(xmatrix<double>& M, double& tol) { //DX20190215 -added tol
    //DX20190215 [OBSOLETE] double tol = _SYM_TOL_;
    uint lead = 1;
    uint rowCount = M.urows;
    uint columnCount = M.ucols;
    xmatrix<double> tmp(rowCount, columnCount);
    for (uint r = 0; r < rowCount; r++) {
      if(columnCount <= lead) {
        break;
      }
      uint i = r;
      while (aurostd::abs(M(i + 1, lead)) < tol) {
        i = i + 1;
        if(rowCount == i) {
          i = r;
          lead = lead + 1;
          if(columnCount == lead) {
            break;
          }
        }
      }
      //Problem with lead +1 here
      if(lead != columnCount) {
        swap_rows(M, i + 1, r + 1);
        //bool consider_half = false;
        double divisor = M(r + 1, lead);
        if(aurostd::abs(M(r + 1, lead)) > tol) {
          for (uint ii = 1; ii <= columnCount; ii++) {
            M(r + 1, ii) = double(M(r + 1, ii) / divisor);
          }
        } else {
          for (uint row = r + 1; row < rowCount; row++) {
            if(aurostd::abs(M(row, lead)) > tol) {
              swap_rows(M, r + 1, row);
              break;
            }
          }
          double divisor = M(r + 1, lead);
          for (uint ii = 1; ii <= columnCount; ii++) {
            M(r + 1, ii) = double(M(r + 1, ii) / divisor);
          }
        }
        for (uint i = 0; i < rowCount; i++) {
          if(i != r) {
            for (uint c = 1; c <= columnCount; c++) {
              tmp(i + 1, c) = M(i + 1, c) - M(i + 1, lead) * M(r + 1, c);
            }
            for (uint c = 1; c <= columnCount; c++) {
              M(i + 1, c) = tmp(i + 1, c);
            }
          }
        }
        lead = lead + 1;
      }
    }
  }
} //namespace SYM

// ******************************************************************************
// find_solution_UD (Find Solution for Underdetermined System)
// ******************************************************************************
namespace SYM {
  bool find_solution_UD(xmatrix<double> M, xvector<double>& SOL, double& tol) { //DX20190215 -added tol
    //M should be in row reduced form, simply check that there are no rows with zero entries and nonzero end value
    //DX20190215 [OBSOLETE] double tol = _SYM_TOL_;
    uint numCols = M.ucols;
    uint numRows = M.urows;
    xvector<double> tmp;
    bool allzero = false;
    for (uint i = 1; i <= numRows; i++) {
      allzero = true;
      for (uint j = 1; j <= numCols - 1; j++) {
        if(aurostd::abs(M(i, j)) > tol) {
          allzero = false;
        }
      }
      if(allzero == true) {
        if(aurostd::abs(M(i, numCols)) > tol && aurostd::abs(M(i, numCols)) < 1 - tol) {
          return false;
        }
      }
    }
    for (uint i = 1; i <= numRows; i++) {
      for (uint j = 1; j <= numCols - 1; j++) {
        if(aurostd::abs(M(i, j) - 1) < tol) {
          tmp(j) = M(i, numCols);
          break;
        }
      }
    }

    SOL = tmp;
    return true;
  }
} //namespace SYM

// ******************************************************************************
// Check Linear System
// ******************************************************************************
// Check for any rows which contain null vectors and a nonzero augmented part
namespace SYM {
  bool checkLinearSystem(vector<xvector<double> >& LHS, vector<double>& RHS, xmatrix<double>& lattice, double& tol) { //DX20190215 -added tol
    //double tol = 5e-3;
    //DX20190215 [OBSOLETE] double tol = _SYM_TOL_;
    bool allzero = true;
    uint numRows = LHS.size();
    uint numCols = 3;

    for (uint i = 0; i < numRows; i++) {
      allzero = true;
      int xyz = 0;
      for (uint j = 1; j <= numCols; j++) {
        if(aurostd::abs(LHS[i][j]) > tol) {
          allzero = false;
          break;
        }
      }
      if(allzero == true) {
        xyz = i % 3 + 1;
        if(aurostd::abs(RHS[i]) > tol) {
          // Check if adding lattice vector components makes the RHS zero.
          if(aurostd::abs(aurostd::abs(RHS[i]) - aurostd::abs(lattice(1)[xyz])) < tol ||
              aurostd::abs(aurostd::abs(RHS[i]) - aurostd::abs(lattice(2)[xyz])) < tol ||
              aurostd::abs(aurostd::abs(RHS[i]) - aurostd::abs(lattice(3)[xyz])) < tol) {
            continue;
          }
          return false;
        }
      }
    }

    // Check if two rows are same, but corresponding shifts (RHS) are different
    for (uint v = 0; v < LHS.size(); v++) {
      int xyz_v = v % 3 + 1;
      for (uint u = v + 1; u < LHS.size(); u++) {
        int xyz_u = u % 3 + 1;
        if(aurostd::abs(LHS[v][1] - LHS[u][1]) < tol && aurostd::abs(LHS[v][2] - LHS[u][2]) < tol && aurostd::abs(LHS[v][3] - LHS[u][3]) < tol) {
          if(aurostd::abs(RHS[v] - RHS[u]) > tol) {
            if(aurostd::abs(aurostd::abs(RHS[v] - RHS[u]) - aurostd::abs(lattice(1)[xyz_v])) < tol ||
                aurostd::abs(aurostd::abs(RHS[v] - RHS[u]) - aurostd::abs(lattice(2)[xyz_v])) < tol ||
                aurostd::abs(aurostd::abs(RHS[v] - RHS[u]) - aurostd::abs(lattice(3)[xyz_v])) < tol ||
                aurostd::abs(aurostd::abs(RHS[v] - RHS[u]) - aurostd::abs(lattice(1)[xyz_u])) < tol ||
                aurostd::abs(aurostd::abs(RHS[v] - RHS[u]) - aurostd::abs(lattice(2)[xyz_u])) < tol ||
                aurostd::abs(aurostd::abs(RHS[v] - RHS[u]) - aurostd::abs(lattice(3)[xyz_u])) < tol) {
              continue;
            }
            return false;
          }
        }
      }
    }
    return true;
  }
} //namespace SYM

// ******************************************************************************
// solve_overdetermined_system (Solve Overdetermined System)
// ******************************************************************************
namespace SYM {
  bool solve_overdetermined_system(vector<xvector<double> >& LHS, vector<double>& RHS, xvector<double>& SOL, xmatrix<double>& lattice, double& min_dist, double& tol) { //DX20190215 - added tol
    //cerr << "SOLVING SYSTEM" << endl;
    //xb();
    //print(LHS);
    //cerr << "-------------------------" << endl;
    //print(RHS);
    //xb();
    //xmatrix<double> coeff_xmat = xvec2xmat(LHS,RHS);
    //xmatrix<double> LHS_xmat = xvec2xmat(LHS);
    //ReducedRowEchelonForm(coeff_xmat,tol);
    //xb();
    //aurostd::cematrix cemat(LHS_xmat);
    //cerr << cemat.InverseMatrix() << endl;;

    xmatrix<double> f2c = trasp(lattice);
    xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = SYM::isLatticeSkewed(lattice, min_dist, tol); //DX20190215 

    vector<xvector<double> > LHS_cart;
    vector<double> RHS_cart;
    for (uint i = 0; i < LHS.size(); i += 3) {
      // Determine Cartesian representation of rotation matrix
      xvector<double> a;
      xvector<double> b;
      xvector<double> c;
      a = LHS[i];
      b = LHS[i + 1];
      c = LHS[i + 2];
      xmatrix<double> lhs_tmp = xvec2xmat(a, b, c);
      xmatrix<double> lhs_cart_mat = f2c * lhs_tmp * inverse(f2c);
      LHS_cart.push_back(lhs_cart_mat(1));
      LHS_cart.push_back(lhs_cart_mat(2));
      LHS_cart.push_back(lhs_cart_mat(3));

      // Determine Cartesian representation of shift
      xvector<double> rhs_tmp;
      rhs_tmp(1) = RHS[i];
      rhs_tmp(2) = RHS[i + 1];
      rhs_tmp(3) = RHS[i + 2];
      xvector<double> rhs_cart_xvec = f2c * rhs_tmp;
      if(skew) {
        xvector<double> origin;
        rhs_cart_xvec = SYM::minimizeDistanceCartesianMethod(origin, rhs_cart_xvec, lattice); //DX20190613
        //DX20190613 [OBSOLETE] SYM::minimizeCartesianDistance(origin, rhs_cart_xvec, rhs_tmp, c2f, f2c, tol); //DX20190215
      } else {
        xvector<double> rhs_frac_xvec = SYM::minimizeDistanceFractionalMethod(rhs_tmp); //DX20190613
        rhs_cart_xvec = f2c * rhs_frac_xvec; //DX20190718 - rhs_frac_xvec not rhs_tmp //DX20190613
        //DX20190613 [OBSOLETE] SYM::PBC(rhs_tmp);
      }
      //DX20190613 [OBSOLETE] rhs_cart_xvec = f2c * rhs_tmp;
      RHS_cart.push_back(rhs_cart_xvec(1));
      RHS_cart.push_back(rhs_cart_xvec(2));
      RHS_cart.push_back(rhs_cart_xvec(3));
    }

    //DX20190215 [OBSOLETE] double tol = _SYM_TOL_;
    xmatrix<double> firstLI;
    int f1, f2, f3;
    bool underdetermined = true;
    //extract three linearly independent vectors
    for (uint i = 0; i < LHS.size(); i++) {
      for (uint j = 1; j < LHS.size(); j++) {
        for (uint k = 2; k < LHS.size(); k++) {
          firstLI = xvec2xmat(LHS[i], LHS[j], LHS[k]);
          if(aurostd::abs(aurostd::determinant(firstLI)) > 1e-10) {
            underdetermined = false;
            f1 = i;
            f2 = j;
            f3 = k;
            i = LHS_cart.size();
            j = LHS_cart.size();
            k = LHS_cart.size();
          }
        }
      }
    }
    if(underdetermined == true) {
      //cerr << "underdetermined" << endl;
      //exists some freedom in selection of origin
      if(!checkLinearSystem(LHS_cart, RHS_cart, lattice,tol)) { //DX20190215
        //cerr << "NULL VECTORS GIVE NONZERO AUGMENTED PART!!!" << endl;
        return false;
      }
      // Need to check different cells (add integers to RHS) since this is a periodic linear system
      for (uint i = 0; i < 3; i++) {
        for (uint j = 0; j < 3; j++) {
          for (uint k = 0; k < 3; k++) {
            vector<double> dS_mod = RHS;
            dS_mod[0] = RHS[0] + double(i);
            dS_mod[1] = RHS[1] + double(j);
            dS_mod[2] = RHS[2] + double(k);
            xmatrix<double> coeff_xmat = xvec2xmat(LHS, dS_mod);
            //cerr << "coeff (before): " << endl;
            //cerr << coeff_xmat << endl;
            ReducedRowEchelonForm(coeff_xmat,tol);
            //cerr << "coeff (after): " << endl;
            //cerr << coeff_xmat << endl;
            if(find_solution_UD(coeff_xmat, SOL, tol)) { //DX20190215
              //cerr << "found sol: " << SOL << endl;
              return true;
            }
          }
        }
      }
      return false;
    }

    xvector<double> dS;
    dS(1) = RHS[f1];
    dS(2) = RHS[f2];
    dS(3) = RHS[f3];

    // Need to check different cells (add integers to RHS) since this is a periodic linear system
    for (uint i = 0; i < 3; i++) {
      for (uint j = 0; j < 3; j++) {
        for (uint k = 0; k < 3; k++) {
          xvector<double> dS_mod;
          dS_mod(1) = dS(1) + double(i);
          dS_mod(2) = dS(2) + double(j);
          dS_mod(3) = dS(3) + double(k);
          //Do you need to put a pseudoinverse here (use SVD)
          xvector<double> origin_shift = aurostd::inverse(firstLI) * dS_mod;
          bool found = true;
          for (uint ii = 0; ii < LHS_cart.size(); ii += 3) {
            xvector<double> coord;
            coord(1) = aurostd::scalar_product(LHS[ii], origin_shift) - RHS[ii]; //DX20200724 - SYM::DotPro to aurostd::scalar_product
            coord(2) = aurostd::scalar_product(LHS[ii + 1], origin_shift) - RHS[ii + 1]; //DX20200724 - SYM::DotPro to aurostd::scalar_product
            coord(3) = aurostd::scalar_product(LHS[ii + 2], origin_shift) - RHS[ii + 2]; //DX20200724 - SYM::DotPro to aurostd::scalar_product
            double min_dist = 1e9; //DX20190613
            if(skew) {
              xvector<double> origin;
              xvector<double> coord_cart = f2c * coord;
              xvector<double> min_coord_cart = SYM::minimizeDistanceCartesianMethod(origin, coord_cart, lattice); //DX20190613
              min_dist = aurostd::modulus(min_coord_cart); //DX20190613
              //DX20190613 [OBSOLETE] if(!SYM::minimizeCartesianDistance(origin, coord_cart, coord, c2f, f2c, tol)) { //DX20190215
              //DX20190613 [OBSOLETE]   found = false;
              //DX20190613 [OBSOLETE]   break;
              //DX20190613 [OBSOLETE] }
            } else {
              xvector<double> min_frac_coord = SYM::minimizeDistanceFractionalMethod(coord); //DX20190613
              min_dist = aurostd::modulus(f2c * min_frac_coord); //DX20190613
              //DX20190613 [OBSOLETE] SYM::PBC(coord);
              //DX20190613 [OBSOLETE] if(aurostd::modulus(f2c * coord) > tol) {
              //DX20190613 [OBSOLETE]   found = false;
              //DX20190613 [OBSOLETE]   break;
              //DX20190613 [OBSOLETE] }
            }
            //DX20190613 - move loop outside if-statements - START
            if(min_dist > tol){
              found = false;
              break;
            }
            //DX20190613 - move loop outside if-statements - END
          }
          if(found == true) {
            //cerr << "FOUND CART" << endl;
            SOL = origin_shift;
            return true;
          }
        }
      }
    }
    return false;
  }
} //namespace SYM

// ******************************************************************************
// get_symmetry_symbols
// ******************************************************************************
//Need effective way to gather all permutations of appropriate combinations (e.g., 2,1,1)
//Function to extract symmetry symbols from wyckoff data (vector<string>)
namespace SYM {
  vector<string> get_symmetry_symbols(string sg) {
    vector<string> out;
    int line = 0;
    //bool add;
    ostringstream tmp;
    tmp.str("");

    for (uint i = 0; i < sg.length(); i++) {
      if(sg[i] == 0x0a) {
        i++;
        line++;
        //add = true;
        //tmp.clear();
      }
      if(line > 1)  // && add == true)
      { //CO20200106 - patching for auto-indenting
        if(sg[i] != '(') {
          tmp << sg[i];
          //cerr << sg[i] << " ";
        } else {
          out.push_back(tmp.str());
          tmp.str("");
          //add = false;
          while (sg[i + 1] != 0x0a) {
            i++;
          }
        }
      }
    }
    return out;
  }
} //namespace SYM

// ******************************************************************************
// get_multiplicities (Obtain Wyckoff Multiplicities from ITC)
// ******************************************************************************
//Function to extract multiplicites from wyckoff data (vector<int>)
namespace SYM {
  vector<int> get_multiplicities(string sg) {
    vector<int> wyckoff_multiplicities;
    int line = 0;
    int mult_factor = 0;
    for (uint i = 0; i < sg.length(); i++) {
      if(sg[i] == 0x0a)
        line++;
      if(line == 2)
        break;
      if(sg[i] == '(')
        mult_factor++;
    }
    line = 0;
    ostringstream temp;
    for (uint i = 0; i < sg.length(); i++) {
      if(sg[i] == 0x0a)
        line++;
      if(line > 1) {
        ostringstream mult;
        temp << sg[i];
        if(sg[i] == 0x0a) {
          temp << temp.str() << endl;
          if(temp.str()[3] != ' ') {
            mult << temp.str()[1] << temp.str()[2] << temp.str()[3] << endl;
            wyckoff_multiplicities.push_back(atoi(mult.str().c_str()));
          } else if(temp.str()[2] != ' ') {
            mult << temp.str()[1] << temp.str()[2] << endl;
            wyckoff_multiplicities.push_back(atoi(mult.str().c_str()));
          } else {
            mult << temp.str()[1] << endl;
            //	  cerr << atoi(mult.str().c_str()) << endl;
            wyckoff_multiplicities.push_back(atoi(mult.str().c_str()));
          }
          temp.str("");
        }
      }
    }
    return wyckoff_multiplicities;
  }
} //namespace SYM

// ******************************************************************************
// get_centering (Obtain Centerings from ITC)
// ******************************************************************************
namespace SYM {
  vector<vector<string> > get_centering(string spaceg) {
    vector<vector<string> > centering;
    vector<string> triplet;
    int i = 0;
    //cerr << spaceg << endl;

    while (spaceg[i] != '(') {  //Get past Space Group #
      //cerr << spaceg[i];
      i++;
    }
    ostringstream os;
    while (spaceg[i] != '\n') {
      //    if(spaceg[i]==' '){
      //      cerr << os.str();
      //      triplet.push_back(os.str());
      //      os.str("");
      //    }
      //cerr << spaceg[i+1] << endl;
      if(spaceg[i] != ',' && spaceg[i] != ')' && spaceg[i] != '(' && spaceg[i] != ' ') {
        if(spaceg[i] == '.') {
          os << spaceg[i];
          i++;
          os << spaceg[i];
          triplet.push_back(os.str());
        } else if(spaceg[i + 1] == '/') {
          os << spaceg[i];
          os << spaceg[i + 1];
          os << spaceg[i + 2];
          i += 2;
          triplet.push_back(os.str());
        } else {
          os << spaceg[i];
          triplet.push_back(os.str());
        }
      }
      os.str("");
      if(triplet.size() == 3) {  //3 Dimensions
        centering.push_back(triplet);
        triplet.clear();
      }
      i++;
    }
    // for(int i=0;i<centering.size();i++){ //OUTPUT
    //   cerr << " " << endl;                       //OUTPUT
    //   for(int j=0;j<3;j++){              //OUTPUT
    //     cerr << centering[i][j];         //OUTPUT
    //   }                                  //OUTPUT
    // }                                     //OUTPUT
    return centering;
  }
} //namespace SYM

// ******************************************************************************
// minimum_enumerated_Wyckoff_letters (Identify minimum enumerated Wyckoff letters) 
// ******************************************************************************
// Given a set of Wyckoff positions, determine the minimum enumerated Wyckoff 
// letter representation that preserves the Wyckoff multiplicity and site 
// symmetry
namespace SYM{
  vector<string> get_minimum_enumerated_Wyckoff_letters(string spacegroupstring, vector<int>& multiplicities, vector<string> site_symmetries){
    vector<string> minimum_enumerated_Wyckoff_letters;
    for (uint i = 0; i < multiplicities.size(); i++) {
      vector<string> site_syms;
      vector<string> letters;
      vector<vector<string> > all_positions;
      // get all Wyckoff positions associated with a given multiplicity and site symmetry 
      //SYM::get_certain_wyckoff_pos(spacegroupstring, multiplicities[i], site_symmetries[i], site_syms, letters, positions);
      SYM::getWyckoffAutomorphismSets(spacegroupstring, multiplicities[i], site_symmetries[i], site_syms, letters, all_positions, false); // false = need to explore lattice axes swaps as well

      vector<int> enumerated_letters;
      int min_enumerated_letter = 1e9;
      string min_letter = "";
      for(uint j=0;j<letters.size();j++){
        int enumerated_letter = SYM::enumerate_wyckoff_letter(letters[j]);
        if(enumerated_letter<min_enumerated_letter){
          bool contains_variable = false;
          if(aurostd::substring2bool(all_positions[j][0], "x") || 
              aurostd::substring2bool(all_positions[j][0], "y") || 
              aurostd::substring2bool(all_positions[j][0], "z")){
            contains_variable = true;
          }
          bool found_minimum_letter = false;
          // check if Wyckoff letter is already used and if it can be used again (i.e., variable coordinate)
          for(uint k=0;k<minimum_enumerated_Wyckoff_letters.size();k++){
            if(letters[j]==minimum_enumerated_Wyckoff_letters[k] && contains_variable){
              min_enumerated_letter=enumerated_letter;
              min_letter=letters[j];
              found_minimum_letter = true;
              break;
            }
            else if(letters[j]==minimum_enumerated_Wyckoff_letters[k] && !contains_variable){
              found_minimum_letter = true;
              break;
            }
          }
          // if letter is not already used, then it is a candidate Wyckoff letter
          if(!found_minimum_letter){
            min_enumerated_letter=enumerated_letter;
            min_letter=letters[j];
          }
        }
      }
      minimum_enumerated_Wyckoff_letters.push_back(min_letter);
    }
    return minimum_enumerated_Wyckoff_letters;
  }
}

//int minimum_enumerated_Wyckoff_sum = SYM::get_minimum_enumerated_Wyckoff_sum(spacegroupstring, wyckoff_mult[w], wyckoff_site_sym[w]); 


// ******************************************************************************
// enumerate_wyckoff_letter (Enumerate the Wyckoff Position Letters)
// ******************************************************************************
namespace SYM {
  int enumerate_wyckoff_letter(string& wyckoff_letter) {
    //First enumerate alphabet
    vector<enum_alpha> alpha_numeric;
    string alphabet = "abcdefghijklmnopqrstuvwxyz";

    //Store alphabet
    for (uint a = 0; a < 26; a++) {
      enum_alpha tmp;
      tmp.letter = alphabet[a];
      tmp.number = a + 1;
      alpha_numeric.push_back(tmp);
    }

    //Store first three greek alphabet
    enum_alpha tmp;
    tmp.letter = "alpha";
    tmp.number = 27;
    alpha_numeric.push_back(tmp);
    tmp.letter = "beta";
    tmp.number = 28;
    alpha_numeric.push_back(tmp);
    tmp.letter = "gamma";
    tmp.number = 29;
    alpha_numeric.push_back(tmp);

    // Determine the corresponding enumeration for the Wyckoff letter input
    int wyckoff_number = 0;
    for (uint a = 0; a < 28; a++) {
      if(wyckoff_letter == alpha_numeric[a].letter) {
        wyckoff_number = alpha_numeric[a].number;
      }
    }
    return wyckoff_number;
  }
} //namespace SYM

// ******************************************************************************
// enumerate_wyckoff_letters (vector) 
// ******************************************************************************
namespace SYM{
  vector<int> enumerate_wyckoff_letters(vector<string>& wyckoff_letters){
    vector<int> wyckoff_numbers;
    for(uint i=0;i<wyckoff_letters.size();i++){
      wyckoff_numbers.push_back(SYM::enumerate_wyckoff_letter(wyckoff_letters[i])); 
    }
    return wyckoff_numbers;
  }
}

// ******************************************************************************
// get_certain_wyckoff_pos (Get certain Wyckoff position based on mult and site sym)
// ******************************************************************************
namespace SYM {
  void get_all_wyckoff_for_site_symmetry(string spaceg, int mult, string site_symmetry, vector<vector<string> >& all_positions) {
    vector<string> sg_lines;
    vector<string> positions;
    stringstream temp;
    int line = 0;
    for (uint i = 0; i < spaceg.length(); i++) {
      if(spaceg[i] == 0x0a)
        line++;
      if(line > 1) {
        temp << spaceg[i];
        if(spaceg[i] == 0x0a) {
          string multiplicity, letter, site_symm, pos;
          temp >> multiplicity;
          temp >> letter;
          temp >> site_symm;
          int tmp_mult = atoi(multiplicity.c_str());
          for(int m=0;m<tmp_mult;m++){
            temp >> pos;
            positions.push_back(pos);
          }
          if(atoi(multiplicity.c_str()) == mult && site_symm == site_symmetry) {
            all_positions.push_back(positions);
          }
          positions.clear();
          temp.str(std::string());
          temp.clear();
        }
      }
    }
  }
} //namespace SYM

// ******************************************************************************
// get_Wyckoff_from_letter
// ******************************************************************************
//DX20191029 END
namespace SYM {
  void get_Wyckoff_from_letter(uint space_group_number, string& space_group_setting, 
      string& Wyckoff_letter, uint& Wyckoff_multiplicity, string& site_symmetry, vector<string>& positions){

    SymmetryInformationITC ITC_sym_info;
    ITC_sym_info.initsgs(space_group_setting);
    string spacegroupstring = ITC_sym_info.gl_sgs[space_group_number - 1];

    get_Wyckoff_from_letter(spacegroupstring, Wyckoff_letter, Wyckoff_multiplicity, site_symmetry, positions);
  }
}
//DX20191029 END

// ---------------------------------------------------------------------------
//DX20190128 START
namespace SYM {
  void get_Wyckoff_from_letter(string& spaceg, string& Wyckoff_letter, 
      uint& Wyckoff_multiplicity, string& site_symmetry, vector<string>& positions) { //DX20190128 - int to uint
    stringstream temp;
    int line = 0;
    for (uint i = 0; i < spaceg.length(); i++) {
      if(spaceg[i] == 0x0a){
        line++;
      }
      if(line > 1) {
        temp << spaceg[i];
        if(spaceg[i] == 0x0a) {
          string multiplicity, letter, site_symm, pos;
          temp >> multiplicity;
          temp >> letter;
          temp >> site_symmetry;
          //DX20190128 [OBSOLETE] - use aurostd - Wyckoff_multiplicity = atoi(multiplicity.c_str());
          Wyckoff_multiplicity = aurostd::string2utype<uint>(multiplicity.c_str());
          if(letter == Wyckoff_letter){
            for(uint m=0;m<Wyckoff_multiplicity;m++){ //DX20190128 - int to uint
              temp >> pos;
              positions.push_back(pos);
            }
            break;
          }
          else {
            temp.str(std::string());
            temp.clear();
            continue;
          }
        }
      }
    }
  }
} //namespace SYM
//DX20190128 END

// ******************************************************************************
// Wyckoff_position_string2vector(string& string_position) 
// ******************************************************************************
// Convert Wyckoff string (from aflow_symmetry_spacegroup_ITC_library) into 
// 
namespace SYM {
  xvector<double> Wyckoff_position_string2xvector(string& string_position){
    xvector<double> vector_position;
    string tmp_string = string_position;
    // clean up string
    aurostd::StringSubst(tmp_string,"(",""); 
    aurostd::StringSubst(tmp_string,")","");
    tmp_string = aurostd::RemoveWhiteSpaces(tmp_string);

    // split string into vector
    vector<string> tokens;
    aurostd::string2tokens(tmp_string,tokens,",");

    // convert vector<string> to xvector<double>
    int coord = 1;
    for(uint i=0;i<tokens.size();i++){
      // check for fraction
      if(aurostd::substring2bool(tokens[i],"/")){
        vector<string> frac_tokens;
        aurostd::string2tokens(tokens[i],frac_tokens,"/");
        if(frac_tokens.size()==2){
          double numerator = aurostd::string2utype<double>(frac_tokens[0]);
          double denominator = aurostd::string2utype<double>(frac_tokens[1]);
          vector_position[coord] = numerator/denominator;
        }
      }
      else {
        vector_position[coord] = aurostd::string2utype<double>(tokens[i]);
      }
      coord++; 
    }

    //stringstream ss_double;
    //string tmp_num;
    //int coord = 1;
    //for (uint i=0; i<string_postition.size(); i++) {
    //if(string_postition[i] != '(' && 
    //string_postition[i] != ',' && 
    //string_postition[i] != ')'){
    //ss_double << string_postition[i];
    //}
    //if(string_postition[i] == ',' || 
    //string_postition[i] == ')'){
    //ss_double >> tmp_num;
    //ss_double.str(std::string());
    //ss_double.clear();
    //// account for Wyckoff position as a fraction
    //if(tmp_num.find("/") != std::string::npos){
    //vector<string> tokens;
    //aurostd::string2tokens(tmp_num,tokens,"/");
    //if(tokens.size()==2){
    //double numerator = aurostd::string2utype<double>(tokens[0]);
    //double denominator = aurostd::string2utype<double>(tokens[1]);
    //vector_position[coord] = numerator/denominator;
    //}
    //}
    //else {
    //vector_position[coord] = aurostd::string2utype<double>(tmp_num);
    //}
    //coord++;
    //}
    //}

    return vector_position;
  }
}

// ******************************************************************************
// get_possible_origin_shifts 
// ******************************************************************************
namespace SYM {
  vector<xvector<double> > get_possible_origin_shifts(string spacegroupstring, int multiplicity, string site_symmetry){
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    vector<xvector<double> > possible_shifts;
    // === Determine possible shifts; get all positions corresponding to the Wyckoff position with the smallest multiplicity in the set and same site symmetry === //
    vector<vector<string> > same_site_symmetry_positions;
    //SYM::get_all_wyckoff_for_site_symmetry(spacegroupstring, multiplicity, site_symmetry, same_site_symmetry_positions);
    vector<string> site_symmetries, letters, positions;
    SYM::getWyckoffAutomorphismSets(spacegroupstring, multiplicity, site_symmetry, site_symmetries, letters, same_site_symmetry_positions, true); // true = check origin shifts only, not axes swaps
    if(LDEBUG) {
      cerr << "SYM::get_possible_origin_shifts: All Wyckoff positions with multiplicty " << multiplicity << " and same site symmetry " << site_symmetry << "." << endl;
      for(uint s=0;s<same_site_symmetry_positions.size();s++){
        print(same_site_symmetry_positions[s]);   
        xb();  
      }
    }
    for (uint s = 0; s < same_site_symmetry_positions.size(); s++) {  //DX20170911 - Check all shifts not just minimum (minimum may not work) 
      for (uint t = 0; t < same_site_symmetry_positions[s].size(); t++) {  //DX20170911 - Check all shifts not just minimum (minimum may not work)
        xvector<double> candidate_shift = SYM::Wyckoff_position_string2xvector(same_site_symmetry_positions[s][t]); 
        bool stored_shift=false; 
        for(uint p=0;p<possible_shifts.size();p++){
          if(aurostd::abs(aurostd::modulus(candidate_shift-possible_shifts[p]))<_ZERO_TOL_){
            stored_shift=true;
            break;
          }
        }
        if(!stored_shift){
          //DX20190905 [OBSOLETE-no more mod_one_xvec] possible_shifts.push_back(SYM::mod_one_xvec(candidate_shift));
          possible_shifts.push_back(BringInCell(candidate_shift)); //DX20190905
          //DX20180228 - Consider difference of two shifts - START
          bool stored_other_shift =false;
          for(uint p=0;p<possible_shifts.size();p++){
            for(uint q=0;q<possible_shifts.size();q++){
              //DX20190905 [OBSOLETE-no more mod_one_xvec] if(aurostd::abs(aurostd::modulus(SYM::mod_one_xvec(candidate_shift-possible_shifts[p]-possible_shifts[q])))<_ZERO_TOL_)
              if(aurostd::abs(aurostd::modulus(BringInCell(candidate_shift-possible_shifts[p]-possible_shifts[q])))<_ZERO_TOL_) //DX20190905
              { //CO20200106 - patching for auto-indenting
                stored_other_shift=true;
                break;
              }
            }
            if(!stored_other_shift){
              //DX20190905 [OBSOLETE-no more mod_one_xvec] possible_shifts.push_back(SYM::mod_one_xvec(candidate_shift-possible_shifts[p]));
              possible_shifts.push_back(BringInCell(candidate_shift-possible_shifts[p])); //DX20190905
            }
          }
          //DX20180228 - Consider difference of two shifts - END
        }
      }
    }
    return possible_shifts;
  }
}

// ******************************************************************************
// get_certain_wyckoff_pos (Get certain Wyckoff position based on mult and site sym)
// ******************************************************************************
namespace SYM {
  void get_certain_wyckoff_pos(string spaceg, int mult, string site_symmetry, vector<string>& site_symmetries, vector<string>& letters, vector<string>& positions) {
    vector<string> sg_lines;
    stringstream temp;
    int line = 0;
    for (uint i = 0; i < spaceg.length(); i++) {
      if(spaceg[i] == 0x0a)
        line++;
      if(line > 1) {
        temp << spaceg[i];
        if(spaceg[i] == 0x0a) {
          string multiplicity, letter, site_symm, first_pos;
          temp >> multiplicity;
          temp >> letter;
          temp >> site_symm;
          temp >> first_pos;
          if(atoi(multiplicity.c_str()) == mult && site_symm == site_symmetry) {
            site_symmetries.push_back(site_symm);
            letters.push_back(letter);
            positions.push_back(first_pos);
          }
          temp.str(std::string());
          temp.clear();
        }
      }
    }
  }
} //namespace SYM

// ******************************************************************************
// getWyckoffAutomorphismSets() //DX20220901
// ******************************************************************************
// Grabs the set of Wyckoff positions that are related by automorphisms of the
// space group. Note: This function uses a trick of arbitrarily ordering the
// Wyckoff site symmetry (the ordering may not be physical). The extracted
// Wyckoff positions may also not be TRULY related by an automorphism; this
// is a greedy algorithm to find the sets that "may" be related by an
// automorphism. Generally, the functions that use this routine should check
// if the sites can actually be mapped onto one another
// The optional boolean "origin_shifts_only" will not consider automorphisms
// related by swapping the lattice labelings
namespace SYM {
  void getWyckoffAutomorphismSets(const string& spaceg, int mult, string site_symmetry, vector<string>& site_symmetries, vector<string>& letters, vector<vector<string> >& all_positions, bool origin_shifts_only) {
   
    vector<string> split_site_symmetry_tmp;
    // store the site symmetry string to match against
    string site_symmetry_base = site_symmetry, site_symmetry_ITC_base = "";
    
    // if exploring all possible automorphisms, including different axes choices,
    // then we need to sort the site symmetry string
    if(!origin_shifts_only) {
      // get sorted Wyckoff site symmetry (note: not necessarily physical, just a speed up trick)
      split_site_symmetry_tmp = SYM::splitSiteSymmetry(site_symmetry);
      std::sort(split_site_symmetry_tmp.begin(),split_site_symmetry_tmp.end());
      site_symmetry_base = aurostd::joinWDelimiter(split_site_symmetry_tmp,"");
    }

    vector<string> all_Wyckoff_strings, Wyckoff_tokens;
    int multiplicity = 0;
    string letter = "", site_symm = "";

    // split up Wyckoff positions
    aurostd::string2tokens(spaceg, all_Wyckoff_strings, "\n");
  
    for(uint i=2; i<all_Wyckoff_strings.size();i++){ //DX 20191107 starting at 2 since the first two lines are header and centering, repectively
      aurostd::string2tokens(all_Wyckoff_strings[i], Wyckoff_tokens, " ");

      // expected sequence: "24 h ..2 (x, y, z) (x, 0, z) ..."
      if(Wyckoff_tokens.size()>3){
        letter = Wyckoff_tokens[1];
        multiplicity = aurostd::string2utype<uint>(Wyckoff_tokens[0]);
        if(multiplicity == mult) {
          site_symm = Wyckoff_tokens[2];
          site_symmetry_ITC_base=site_symm; 

          // if exploring all possible automorphisms, including different axes choices,
          // then we need to sort the site symmetry string
          if(!origin_shifts_only) {
            // sort ITC position
            split_site_symmetry_tmp = SYM::splitSiteSymmetry(site_symm);
            std::sort(split_site_symmetry_tmp.begin(),split_site_symmetry_tmp.end());
            site_symmetry_ITC_base = aurostd::joinWDelimiter(split_site_symmetry_tmp,"");
          }

          if(site_symmetry_ITC_base == site_symmetry_base) {
            site_symmetries.push_back(site_symm);
            letters.push_back(letter);
            vector<string> positions;
            for(uint p=3;p<Wyckoff_tokens.size();p++) {
              positions.push_back(Wyckoff_tokens[p]);
            }
            all_positions.push_back(positions);
          }
        }
      }
    }
  }
} //namespace SYM

// ******************************************************************************
// get_wyckoff_equations (Obtain Wyckoff Equationss from ITC) //DX20170831
// ******************************************************************************
//Function to get wyckoff data in individual lines (vector<string>). The string input is the space group data contained in the wyckoff library. Mult specifies the degeneracy of interest. For example, if you are looking for all wyckoff data with multiplicity 2 in a certain space group.
namespace SYM {
  vector<string> get_wyckoff_equation(string spaceg, int mult) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    vector<int> mult_vec = get_multiplicities(spaceg);
    //Error if mult is not contained in mult_vec (i.e., a wyckoff position with multiplicity mult does not exist for the space group spaceg)
    if(!aurostd::WithinList(mult_vec, mult)) { //DX20210422 - SYM::invec() -> aurostd::WithinList()
      if(LDEBUG) { cerr << "SYM::get_wyckoff_equation: WARNING: no wyckoff position with multiplicity " << mult << "." << endl; }
      vector<string> none;
      return none;
    }
    vector<int> index;
    for (uint i = 0; i < mult_vec.size(); i++) {
      //if(mult_vec[i]*get_centering(spaceg).size()==mult)
      if(mult_vec[i] == mult)
      { //CO20200106 - patching for auto-indenting
        //cerr <<"match: " << i << endl;
        index.push_back(i);
      }
    }

    vector<string> sg_lines;
    ostringstream temp;
    int line = 0;
    for (uint i = 0; i < spaceg.length(); i++) {
      if(spaceg[i] == 0x0a)
        line++;
      if(line > 1) {
        temp << spaceg[i];
        if(spaceg[i] == 0x0a) {
          sg_lines.push_back(temp.str());
          temp.str("");
        }
      }
    }
    vector<string> equations;

    uint start = 0; //DX20180514 - added initialization

    //for(int l=0;l<sg_lines.size();l++)
    for (uint l = index[0]; l < index[0] + index.size(); l++)
    { //CO20200106 - patching for auto-indenting
      for (uint i = 0; i < sg_lines[l].length(); i++) {
        //cerr << sg_lines[l][i] << " ";
        if(sg_lines[l][i] == '(') {
          start = i;
          break;
        }
      }
      ostringstream temp_eqn;
      for (uint i = start; i < sg_lines[l].length(); i++) {
        if(sg_lines[l][i] != '(' && sg_lines[l][i] != ')' && sg_lines[l][i] != ' ') {
          temp_eqn << sg_lines[l][i];
        }
        if(sg_lines[l][i] == ')') {
          equations.push_back(temp_eqn.str());
          temp_eqn.str("");
        }
      }
    }
    return equations;
  }
} //namespace SYM

//DX20190128 START
// ******************************************************************************
// get_wyckoff_pos (Obtain Wyckoff Positions from ITC) (Overloaded)
// ******************************************************************************
//Function to get wyckoff data in individual lines (vector<string>). The string input is the space group data contained in the wyckoff library. Mult specifies the degeneracy of interest. For example, if you are looking for all wyckoff data with multiplicity 2 in a certain space group.
namespace SYM {
  vector<vector<string> > get_wyckoff_pos(string spaceg, int Wyckoff_multiplicity, string Wyckoff_letter) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    vector<int> mult_vec = get_multiplicities(spaceg);
    //Error if mult is not contained in mult_vec (i.e., a wyckoff position with multiplicity mult does not exist for the space group spaceg)
    if(!aurostd::WithinList(mult_vec, Wyckoff_multiplicity)) { //DX20210422 - SYM::invec() -> aurostd::WithinList()
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << ": WARNING: no wyckoff position with multiplicity " << Wyckoff_multiplicity << "." << endl; }
      vector<vector<string> > none;
      return none;
    }

    //DX20191030 - use tokens instead of stringstream assignment - START
    vector<string> all_positions;
    vector<string> all_Wyckoff_strings, Wyckoff_tokens; 
    string multiplicity, letter, site_symmetry, position;

    // split up Wyckoff positions
    aurostd::string2tokens(spaceg, all_Wyckoff_strings, "\n");

    for(uint i=2; i<all_Wyckoff_strings.size();i++){ //DX20191107 starting at 2 since the first two lines are header and centering, repectively
      aurostd::string2tokens(all_Wyckoff_strings[i], Wyckoff_tokens, " ");

      // expected sequence: "24 h ..2 (x, y, z) (x, 0, z) ..."
      if(Wyckoff_tokens.size()>3){
        letter = Wyckoff_tokens[1];
        if(letter == Wyckoff_letter){
          multiplicity = Wyckoff_tokens[0]; site_symmetry = Wyckoff_tokens[2];
          for(uint t=3;t<Wyckoff_tokens.size();t++){ //DX20191030 - need loop to get all positions
            string position = aurostd::RemoveWhiteSpaces(Wyckoff_tokens[t]); position = aurostd::RemoveCharacter(position,'\n');
            if(position.size()!=0){
              all_positions.push_back(position);
            }
          }
          break; //DX20200423 - moved break outside of loop
        }
      }
    }
    //DX20191030 - remove stringstream assignment (used to be here)
    //DX20191030 - use tokens instead of stringstream assignment - END

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << ": Wyckoff positions without centering(s): " << aurostd::joinWDelimiter(all_positions," ") << endl;}

    vector<vector<string> > non_shifted_Wyckoff_positions;
    vector<string> tokens;
    //format positions
    for(uint i=0;i<all_positions.size();i++){
      aurostd::string2tokens(all_positions[i],tokens,",");
      vector<string> vec_pos;
      for(uint t=0;t<tokens.size();t++){
        string pos = aurostd::StringSubst(tokens[t],"(",""); pos = aurostd::StringSubst(pos,")","");
        vec_pos.push_back(pos);
      }
      non_shifted_Wyckoff_positions.push_back(vec_pos);
    }

    vector<vector<string> > all_Wyckoff_positions;

    stringstream tmp_Wyckoff_equation;
    for(uint i=0;i<get_centering(spaceg).size();i++) {
      for(uint j=0;j<non_shifted_Wyckoff_positions.size();j++){
        vector<string> vec_pos;
        for(uint k=0;k<non_shifted_Wyckoff_positions[j].size();k++){
          tmp_Wyckoff_equation << non_shifted_Wyckoff_positions[j][k];
          tmp_Wyckoff_equation << "+" << get_centering(spaceg)[i][k];
          vec_pos.push_back(tmp_Wyckoff_equation.str());
          tmp_Wyckoff_equation.str(std::string());
          tmp_Wyckoff_equation.clear();
        }
        all_Wyckoff_positions.push_back(vec_pos);
      }
    }

    if(LDEBUG) {
      vector<string> tmp; 
      for(uint i=0;i<all_Wyckoff_positions.size();i++){tmp.push_back("("+aurostd::joinWDelimiter(all_Wyckoff_positions[i],",")+")");}
      cerr << __AFLOW_FUNC__ << ": All Wyckoff positions: " << aurostd::joinWDelimiter(tmp," ") << endl;
    }

    return all_Wyckoff_positions;
  }
}
//DX20190128 END

/////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS: get Wyckoff position information (DX20191030)
//
// The following inputs are required to get the Wyckoff position information (minimum)
//   1) space group number
//   2) space group setting
//   3) Wyckoff letter designation
// Other information - Wyckoff multiplicity and site symmetry - are not unique 
// to a given position
//
// Functions below get the equation, multiplicity, site symmetry, or all information

// ******************************************************************************
// getWyckoffEquations 
// ******************************************************************************
namespace SYM {
  // ---------------------------------------------------------------------------
  // input: space group number + setting
  vector<vector<string> > getWyckoffEquations(const uint space_group_number, const string& space_group_setting, const string& Wyckoff_letter){

    // get Wyckoff information based on space group and setting
    SymmetryInformationITC ITC_sym_info;
    ITC_sym_info.initsgs(space_group_setting);
    string Wyckoff_string = ITC_sym_info.gl_sgs[space_group_number - 1];

    return getWyckoffEquations(Wyckoff_string, Wyckoff_letter);
  }

  // ---------------------------------------------------------------------------
  // input: Wyckoff string 
  vector<vector<string> > getWyckoffEquations(const string& Wyckoff_string, const string& Wyckoff_letter){

    uint Wyckoff_multiplicity = 0;
    string site_symmetry = "";
    vector<vector<string> > positions;

    getWyckoffInformation(Wyckoff_string, Wyckoff_letter, Wyckoff_multiplicity, site_symmetry, positions);

    return positions;
  }
}

// ******************************************************************************
// getWyckoffMultiplicity 
// ******************************************************************************
namespace SYM {
  // ---------------------------------------------------------------------------
  // input: space group number + setting
  uint getWyckoffMultiplicity(const uint space_group_number, const string& space_group_setting, const string& Wyckoff_letter){

    // get Wyckoff information based on space group and setting
    SymmetryInformationITC ITC_sym_info;
    ITC_sym_info.initsgs(space_group_setting);
    string Wyckoff_string = ITC_sym_info.gl_sgs[space_group_number - 1];

    return getWyckoffMultiplicity(Wyckoff_string, Wyckoff_letter);
  }

  // ---------------------------------------------------------------------------
  // input: Wyckoff string 
  uint getWyckoffMultiplicity(const string& Wyckoff_string, const string& Wyckoff_letter){

    uint Wyckoff_multiplicity = 0;
    string site_symmetry = "";
    vector<vector<string> > positions;

    getWyckoffInformation(Wyckoff_string, Wyckoff_letter, Wyckoff_multiplicity, site_symmetry, positions);

    return Wyckoff_multiplicity;
  }
}

// ******************************************************************************
// getWyckoffSiteSymmetry
// ******************************************************************************
namespace SYM {
  // ---------------------------------------------------------------------------
  // input: space group number + setting
  string getWyckoffSiteSymmetry(const uint space_group_number, const string& space_group_setting, const string& Wyckoff_letter){

    // get Wyckoff information based on space group and setting
    SymmetryInformationITC ITC_sym_info;
    ITC_sym_info.initsgs(space_group_setting);
    string Wyckoff_string = ITC_sym_info.gl_sgs[space_group_number - 1];

    return getWyckoffSiteSymmetry(Wyckoff_string, Wyckoff_letter);
  }

  // ---------------------------------------------------------------------------
  // input: Wyckoff string 
  string getWyckoffSiteSymmetry(const string& Wyckoff_string, const string& Wyckoff_letter){

    uint Wyckoff_multiplicity = 0;
    string site_symmetry = "";
    vector<vector<string> > positions;

    getWyckoffInformation(Wyckoff_string, Wyckoff_letter, Wyckoff_multiplicity, site_symmetry, positions);

    return site_symmetry;
  }
}

// ******************************************************************************
// getWyckoffInformation 
// ******************************************************************************
namespace SYM {
  // ---------------------------------------------------------------------------
  // input: space group number + setting
  void getWyckoffInformation(const uint space_group_number, const string& space_group_setting, const string& Wyckoff_letter,
      uint& Wyckoff_multiplicity, string& site_symmetry, vector<vector<string> >& all_positions){

    // get Wyckoff information based on space group and setting
    SymmetryInformationITC ITC_sym_info;
    ITC_sym_info.initsgs(space_group_setting);
    string Wyckoff_string = ITC_sym_info.gl_sgs[space_group_number - 1];

    return getWyckoffInformation(Wyckoff_string, Wyckoff_letter, Wyckoff_multiplicity, site_symmetry, all_positions);
  }

  // ---------------------------------------------------------------------------
  // input: space group number + setting
  void getWyckoffInformation(const string& Wyckoff_string, const string& Wyckoff_letter,
      uint& Wyckoff_multiplicity, string& site_symmetry, vector<vector<string> >& all_positions){

    bool LDEBUG = (FALSE || XHOST.DEBUG);
    bool reduce = true; // DEBUGGING variable: simplify/reduce Wyckoff coordinates (e.g., 0.25+0.5 -> 0.75)
    vector<string> split_Wyckoff_strings, Wyckoff_tokens, positions;

    // split up Wyckoff positions
    aurostd::string2tokens(Wyckoff_string, split_Wyckoff_strings, "\n");

    for(uint i=2; i<split_Wyckoff_strings.size();i++){ //DX20191107 starting at 2 since the first two lines are header and centering, repectively
      aurostd::string2tokens(split_Wyckoff_strings[i], Wyckoff_tokens, " ");

      // expected sequence: "24 h ..2 (x, y, z) (x, 0, z) ..."
      if(Wyckoff_tokens.size()>3){
        string letter = Wyckoff_tokens[1];
        if(letter == Wyckoff_letter){
          Wyckoff_multiplicity = aurostd::string2utype<uint>(Wyckoff_tokens[0]); site_symmetry = Wyckoff_tokens[2];
          for(uint t=3;t<Wyckoff_tokens.size();t++){ //DX20191030 - need loop to get all positions
            string position = aurostd::RemoveWhiteSpaces(Wyckoff_tokens[t]); position = aurostd::RemoveCharacter(position,'\n');
            if(position.size()!=0){
              positions.push_back(position);
            }
          }
          break; //DX20200423 - moved break outside of loop
        }
      }
    }


    if(LDEBUG) {cerr << __AFLOW_FUNC__ << ": Wyckoff positions without centering(s): " << aurostd::joinWDelimiter(positions," ") << endl;}

    vector<vector<string> > non_shifted_Wyckoff_positions;
    vector<string> tokens;
    //format positions
    for(uint i=0;i<positions.size();i++){
      aurostd::string2tokens(positions[i],tokens,",");
      vector<string> vec_pos;
      for(uint t=0;t<tokens.size();t++){
        string pos = aurostd::StringSubst(tokens[t],"(",""); pos = aurostd::StringSubst(pos,")","");
        // ---------------------------------------------------------------------------
        // split into equation entities (number, variable)
        if(reduce){
          vector<sdouble> sd_coordinate = simplify(pos);
          vec_pos.push_back(SYM::formatWyckoffPosition(sd_coordinate));
        }
        else{
          vec_pos.push_back(pos);
        }
      }
      non_shifted_Wyckoff_positions.push_back(vec_pos);
    }

    stringstream tmp_Wyckoff_equation;
    tmp_Wyckoff_equation.str(std::string());
    tmp_Wyckoff_equation.clear();
    for(uint i=0;i<get_centering(Wyckoff_string).size();i++) {
      for(uint j=0;j<non_shifted_Wyckoff_positions.size();j++){
        vector<string> vec_pos;
        for(uint k=0;k<non_shifted_Wyckoff_positions[j].size();k++){
          tmp_Wyckoff_equation << non_shifted_Wyckoff_positions[j][k];
          tmp_Wyckoff_equation << "+" << get_centering(Wyckoff_string)[i][k];
          // ---------------------------------------------------------------------------
          // split into equation entities (number, variable) 
          if(reduce){
            vector<sdouble> sd_coordinate = simplify(tmp_Wyckoff_equation.str());
            vec_pos.push_back(SYM::formatWyckoffPosition(sd_coordinate));
          }
          else{
            vec_pos.push_back(tmp_Wyckoff_equation.str());
          }
          tmp_Wyckoff_equation.str(std::string());
          tmp_Wyckoff_equation.clear();
        }
        all_positions.push_back(vec_pos);
      }
    }

    if(LDEBUG) {
      vector<string> tmp_vposition; 
      for(uint i=0;i<all_positions.size();i++){tmp_vposition.push_back("("+aurostd::joinWDelimiter(all_positions[i],",")+")");}
      cerr << __AFLOW_FUNC__ << ": All Wyckoff positions: " << aurostd::joinWDelimiter(tmp_vposition," ") << endl;
    }

  }
}
//DX20190128 END

// ******************************************************************************
// get_wyckoff_pos (Obtain Wyckoff Positions from ITC) (Overloaded)
// ******************************************************************************
//Function to get wyckoff data in individual lines (vector<string>). The string input is the space group data contained in the wyckoff library. Mult specifies the degeneracy of interest. For example, if you are looking for all wyckoff data with multiplicity 2 in a certain space group.
namespace SYM {
  vector<vector<vector<string> > > get_wyckoff_pos(string spaceg, int mult) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    vector<int> mult_vec = get_multiplicities(spaceg);
    //Error if mult is not contained in mult_vec (i.e., a wyckoff position with multiplicity mult does not exist for the space group spaceg)
    if(!aurostd::WithinList(mult_vec, mult)) { //DX20210422 - SYM::invec() -> aurostd::WithinList()
      if(LDEBUG) { cerr << "SYM::get_wyckoff_pos: WARNING: no wyckoff position with multiplicity " << mult << "." << endl; }
      vector<vector<vector<string> > > none;
      return none;
    }
    vector<int> index;
    for (uint i = 0; i < mult_vec.size(); i++) {
      //if(mult_vec[i]*get_centering(spaceg).size()==mult)
      if(mult_vec[i] == mult)
      { //CO20200106 - patching for auto-indenting
        //cerr <<"match: " << i << endl;
        index.push_back(i);
      }
    }

    vector<string> sg_lines;
    ostringstream temp;
    int line = 0;
    for (uint i = 0; i < spaceg.length(); i++) {
      if(spaceg[i] == 0x0a)
        line++;
      if(line > 1) {
        temp << spaceg[i];
        if(spaceg[i] == 0x0a) {
          sg_lines.push_back(temp.str());
          temp.str("");
        }
      }
    }
    vector<string> singlets;
    vector<vector<string> > tuplets;
    //vector<vector<vector<string> > > wyckoff_positions; //DX20170831
    vector<vector<vector<string> > > wyckoff_set; //DX20170831
    vector<vector<vector<vector<string> > > > all_wyckoff_sets; //DX20170831

    uint start = 0; //DX20180514 - added initialization

    //for(int l=0;l<sg_lines.size();l++)
    for (uint l = index[0]; l < index[0] + index.size(); l++)
    { //CO20200106 - patching for auto-indenting
      for (uint i = 0; i < sg_lines[l].length(); i++) {
        //cerr << sg_lines[l][i] << " ";
        if(sg_lines[l][i] == '(') {
          start = i;
          break;
        }
      }
      ostringstream temp1;
      for (uint i = start; i < sg_lines[l].length(); i++) {
        //if(sg_lines[1][i]=='('){
        //temp1;
        //}
        if(sg_lines[l][i] != '(' && sg_lines[l][i] != ')' && sg_lines[l][i] != ',' && sg_lines[l][i] != ' ') {
          temp1 << sg_lines[l][i];
          singlets.push_back(temp1.str());
          temp1.str("");
        }
        if(sg_lines[l][i] == ',' || sg_lines[l][i] == ')') {
          //cerr << temp1.str() << endl;
          temp1.str("");
          tuplets.push_back(singlets);
          singlets.clear();
        }
        if(sg_lines[l][i] == ')') {
          temp1.str("");
          wyckoff_set.push_back(tuplets); //DX20170831 - wyckoff_positions to wyckoff_set
          tuplets.clear();
        }
      }
      all_wyckoff_sets.push_back(wyckoff_set); //DX20170831
      wyckoff_set.clear(); //DX20170831
    }

    //DISPLAY CENTRING
    //for(uint ii=0;ii<get_centering(spaceg).size();ii++){
    //  for(uint j=0;j<get_centering(spaceg)[ii].size();j++){
    //    cerr << get_centering(spaceg)[ii][j] << " ";
    //  }
    //  cerr << endl;
    //}

    //Wyckoff_positions[tuplets][singlets][singlet characters]
    //First loop is to account for centering operations.

    //Something going on here /db
    //OUTPUT
    vector<string> singletvec;
    vector<vector<string> > tupletvec;
    vector<vector<vector<string> > > outvec;
    ostringstream oss;
    int reg = 0;   //Counter for tuplets
    int breg = 0;  //Counter for multiplicity blocks

    //DX20170831 - Altered scheme to store unshifted positions first, then the shifted positions -START
    for (uint w = 0; w < all_wyckoff_sets.size(); w++){
      for (uint ii = 0; ii < get_centering(spaceg).size(); ii++) {
        for (uint k = 0; k < all_wyckoff_sets[w].size(); k++) {
          for (uint j = 0; j < all_wyckoff_sets[w][k].size(); j++)  //for(int j=0;j<3;j++)
          { //CO20200106 - patching for auto-indenting
            for (uint i = 0; i < all_wyckoff_sets[w][k][j].size(); i++) {
              if(all_wyckoff_sets[w][k][j][i] != "\n")
                oss << all_wyckoff_sets[w][k][j][i];
            }
            oss << "+" << get_centering(spaceg)[ii][j];
            singletvec.push_back(oss.str());
            oss.str("");
            reg++;
            if(reg % 3 == 0) {  //New line every tuplet
              //oss << endl;
              tupletvec.push_back(singletvec);
              //cerr << singletvec.size() << endl;
              singletvec.clear();
              breg++;
            }
          }
          if(breg % mult == 0) {  //New block for every wyckoff position
            outvec.push_back(tupletvec);
            tupletvec.clear();
            //oss << endl;
          }
        }
      }
    }
    //DX20170831 - Altered scheme to store unshifted positions first, then the shifted positions -END
    //cerr << "outvec.size(): " << outvec.size() << endl;
    //for (uint l=0;l<outvec.size(); l++){
    //	cerr << "outvec[l].size(): " << outvec[l].size() << endl;
    //  for (uint m=0;m<outvec[l].size(); m++){
    //    cerr << "outvec[l][m].size(): " << outvec[l][m].size() << endl;
    //    for (uint n=0;n<outvec[l][m].size(); n++){
    //      cerr << "outvec: (l,m,n)" << l << m << n << " " << outvec[l][m][n] << endl;
    //    }
    //  }
    //}
    //cerr << oss.str() << endl;
    //cerr << outvec[0][1][0] << endl;
    //cerr << outvec[0][1][1] << endl;
    //cerr << outvec[0][1][2] << endl;
    return outvec;
    //return wyckoff_positions;
  }
} //namespace SYM

// ******************************************************************************
// get_wyckoff_pos (Obtain Wyckoff Positions from ITC) (Overloaded)
// ******************************************************************************
namespace SYM {
  vector<vector<vector<string> > > get_wyckoff_pos(string spaceg, int mult, bool getcentering) {
    //THE MULT SHOULD BE THE FULL MULTIPLICITY (INCLUDING THE CENTERING OPERATIONS). THIS WILL BE MODIFIED AUTOMATICALLY IF GETCENTERING == FALSE.
    vector<int> mult_vec = get_multiplicities(spaceg);
    //Error if mult is not contained in mult_vec (i.e., a wyckoff position with multiplicity mult does not exist for the space group spaceg)
    if(!aurostd::WithinList(mult_vec, mult)) { //DX20210422 - SYM::invec() -> aurostd::WithinList()
      cerr << "SYM::get_wyckoff_pos: WARNING: no wyckoff position with multiplicity " << mult << "." << endl;
      vector<vector<vector<string> > > none;
      return none;
      ;
    }
    //modify mult and mult_vec if getcentering == false.
    if(getcentering == false) {
      for (uint i = 0; i < mult_vec.size(); i++) {
        mult_vec[i] = mult_vec[i] / get_centering(spaceg).size();
      }
      mult = mult / get_centering(spaceg).size();
    }

    vector<int> index;
    for (uint i = 0; i < mult_vec.size(); i++) {
      if(mult_vec[i] == mult) {
        //cerr <<"match: " << i << endl;
        index.push_back(i);
      }
    }

    vector<string> sg_lines;
    ostringstream temp;
    int line = 0;
    for (uint i = 0; i < spaceg.length(); i++) {
      if(spaceg[i] == 0x0a)
        line++;
      if(line > 1) {
        temp << spaceg[i];
        if(spaceg[i] == 0x0a) {
          sg_lines.push_back(temp.str());
          temp.str("");
        }
      }
    }
    vector<string> singlets;
    vector<vector<string> > tuplets;
    vector<vector<vector<string> > > wyckoff_positions;

    uint start = 0; //DX20180514 - added initialization

    //for(int l=0;l<sg_lines.size();l++){ //[CO20200106 - close bracket for indenting]}
    for (uint l = index[0]; l < index[0] + index.size(); l++) {
      for (uint i = 0; i < sg_lines[l].length(); i++) {
        if(sg_lines[l][i] == '(') {
          start = i;
          break;
        }
      }
      ostringstream temp1;
      for (uint i = start; i < sg_lines[l].length(); i++) {
        //if(sg_lines[1][i]=='('){
        //temp1;
        //}
        if(sg_lines[l][i] != '(' && sg_lines[l][i] != ')' && sg_lines[l][i] != ',' && sg_lines[l][i] != ' ') {
          temp1 << sg_lines[l][i];
          singlets.push_back(temp1.str());
          temp1.str("");
        }
        if(sg_lines[l][i] == ',' || sg_lines[l][i] == ')') {
          //cerr << temp1.str() << endl;
          temp1.str("");
          tuplets.push_back(singlets);
          singlets.clear();
        }
        if(sg_lines[l][i] == ')') {
          temp1.str("");
          wyckoff_positions.push_back(tuplets);
          tuplets.clear();
        }
      }
    }

    //DISPLAY CENTRING
    //  for(uint ii=0;ii<get_centering(spaceg).size();ii++){
    //    for(uint j=0;j<get_centering(spaceg)[ii].size();j++){
    //      cerr << get_centering(spaceg)[ii][j] << " ";
    //    }
    //    cerr << endl;
    //  }

    //Wyckoff_positions[tuplets][singlets][singlet characters]
    //First loop is to account for centering operations.

    //OUTPUT

    //  uint C=1;
    //  if(getcentering == true){
    //    C = get_centering(spaceg).size();
    //  }

    vector<string> singletvec;
    vector<vector<string> > tupletvec;
    vector<vector<vector<string> > > outvec;
    ostringstream oss;
    int reg = 0;   //Counter for tuplets
    int breg = 0;  //Counter for multiplicity blocks
    for (uint k = 0; k < wyckoff_positions.size(); k++) {
      //for(uint ii=0;ii<get_centering(nocent).size();ii++)
      //for(uint ii=0;ii<C;ii++)
      //{ //CO20200106 - patching for auto-indenting
      for (uint j = 0; j < wyckoff_positions[k].size(); j++)  //for(int j=0;j<3;j++)
      { //CO20200106 - patching for auto-indenting
        for (uint i = 0; i < wyckoff_positions[k][j].size(); i++) {
          if(wyckoff_positions[k][j][i] != "\n")
            oss << wyckoff_positions[k][j][i];
        }
        //oss << "+" << get_centering(spaceg)[ii][j];
        singletvec.push_back(oss.str());
        oss.str("");
        reg++;

        if(reg % 3 == 0) {  //New line every tuplet
          //oss << endl;
          tupletvec.push_back(singletvec);
          //cerr << singletvec.size() << endl;
          singletvec.clear();
          breg++;
        }
      }
      if(breg % mult == 0) {  //New block for every wyckoff position
        //cerr << "breg: " << breg << endl;
        outvec.push_back(tupletvec);
        tupletvec.clear();
        //oss << endl;
      }
      //}
    }
    //cerr << oss.str() << endl;
    //cerr << outvec[0][1][0] << endl;
    //cerr << outvec[0][1][1] << endl;
    //cerr << outvec[0][1][2] << endl;
    return outvec;
    //return wyckoff_positions;
  }
} //namespace SYM

// ******************************************************************************
// retunrGeneralWyckoffPosition //DX20170831
// ******************************************************************************
namespace SYM {
  void getGeneralWyckoffMultiplicityAndPosition(uint space_group_number, string& space_group_setting, int& general_wyckoff_multiplicity, vector<string>& general_wyckoff_position){
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    SymmetryInformationITC ITC_sym_info; //DX20190215
    ITC_sym_info.initsgs(space_group_setting); //DX20190215

    //DX20190215 [OBSOLETE] using SYM::gl_sgs;
    string spacegroupstring = ITC_sym_info.gl_sgs[space_group_number - 1];
    vector<int> wyckoff_multiplicities = SYM::get_multiplicities(spacegroupstring);
    general_wyckoff_multiplicity = wyckoff_multiplicities[1];
    general_wyckoff_position = findGeneralWyckoffPosition(spacegroupstring, wyckoff_multiplicities[1]);
    if(LDEBUG) {
      cerr << "SYM::getGeneralWyckoffMultiplicityAndPosition: General Wyckoff position" << endl;
      print(general_wyckoff_position);
    }
  }
}

// ******************************************************************************
// findGeneralWyckoffPosition //DX20170831
// ******************************************************************************
namespace SYM {
  vector<string> findGeneralWyckoffPosition(string& spacegroupstring, int& general_wyckoff_multiplicity){
    vector<string> general_positions;
    vector<vector<vector<string> > > testvvvstring = get_wyckoff_pos(spacegroupstring, general_wyckoff_multiplicity);

    // DEBUG for(uint i=0;i<testvvvstring.size();i++){
    // DEBUG  for(uint j=0;j<testvvvstring[i].size();j++){
    // DEBUG   ::print(testvvvstring[i][j]);
    // DEBUG  }
    //}
    vector<vector<string> > centering = get_centering(spacegroupstring);
    vector<vector<vector<vector<sdouble> > > > testvvvsd = convert_wyckoff_pos_sd(testvvvstring);
    for (uint j = 0; j < testvvvsd.size(); j++) {
      for (uint m = 0; m < testvvvsd[j].size(); m++) {
        string equation = "";
        vector<string> eqn;
        for (uint k = 0; k < 3; k++) { 
          eqn.push_back(SYM::formatWyckoffPosition(testvvvsd[j][m][k])); //DX20190723 - condensed lines below (in previous version) into function
        }
        equation = aurostd::joinWDelimiter(eqn,",");
        general_positions.push_back(equation);
      }
    }
    return general_positions;
  }
} //namespace SYM

//DX20191029 - START
// ******************************************************************************
// wyckoffsite_ITC::getWyckoffFromLetter()
// ******************************************************************************
void wyckoffsite_ITC::getWyckoffFromLetter(uint space_group_number, const string& Wyckoff_letter, int setting){

  SymmetryInformationITC ITC_sym_info;
  string space_group_setting=aurostd::utype2string<int>(setting);
  ITC_sym_info.initsgs(space_group_setting);
  string space_group_string = ITC_sym_info.gl_sgs[space_group_number - 1];

  getWyckoffFromLetter(space_group_string, Wyckoff_letter);
}

void wyckoffsite_ITC::getWyckoffFromLetter(const string& space_group_string, const string& Wyckoff_letter){

  bool LDEBUG = (FALSE || XHOST.DEBUG);
  stringstream message;

  //DX 20191030 - use tokens instead of stringstream assignment - START
  vector<string> all_positions;
  vector<string> all_Wyckoff_strings, Wyckoff_tokens;
  vector<string> components;
  vector<vector<string> > equations_no_centering;
  //string multiplicity, letter, site_symmetry, position;

  // split up Wyckoff positions
  aurostd::string2tokens(space_group_string, all_Wyckoff_strings, "\n");

  for(uint i=2; i<all_Wyckoff_strings.size();i++){ //DX 20191107 starting at 2 since the first two lines are header and centering, repectively
    aurostd::string2tokens(all_Wyckoff_strings[i], Wyckoff_tokens, " ");

    // expected sequence: "24 h ..2 (x, y, z) (x, 0, z) ..."
    if(Wyckoff_tokens.size()>3){
      letter = Wyckoff_tokens[1];
      if(letter == Wyckoff_letter){
        //if(isdigit(Wyckoff_tokens[0])){
        multiplicity = aurostd::string2utype<uint>(Wyckoff_tokens[0]);
        //}
        //else{
        //  message << "The first position of the Wyckoff string is not a number (multiplicity): " << all_Wyckoff_strings[i];
        //  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
        //}
        site_symmetry = Wyckoff_tokens[2];
        for(uint t=3;t<Wyckoff_tokens.size();t++){ //DX 20191030 - need loop to get all positions
          string position = aurostd::RemoveWhiteSpaces(Wyckoff_tokens[t]); position = aurostd::RemoveCharacter(position,'\n');
          if(position.size()!=0){
            components.clear();
            aurostd::string2tokens(position,components,",");
            vector<string> vec_pos;
            for(uint t=0;t<components.size();t++){
              string component = aurostd::StringSubst(components[t],"(",""); component = aurostd::StringSubst(component,")","");
              vec_pos.push_back(component);
            }
            equations_no_centering.push_back(vec_pos);
          }
        }
        break; //DX20200423 - moved break outside of loop
      }
    }
  }
  //DX 20191030 - remove stringstream assignment (used to be here)
  //DX 20191030 - use tokens instead of stringstream assignment - END

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << ": Wyckoff positions without centering(s): ";
    for(uint i=0;i<equations_no_centering.size();i++){
      cerr << aurostd::joinWDelimiter(equations_no_centering[i],",") << endl;
    }
  }

  stringstream tmp_Wyckoff_equation;
  for(uint i=0;i<SYM::get_centering(space_group_string).size();i++) {
    for(uint j=0;j<equations_no_centering.size();j++){
      vector<string> vec_pos;
      for(uint k=0;k<equations_no_centering[j].size();k++){
        tmp_Wyckoff_equation << equations_no_centering[j][k];
        tmp_Wyckoff_equation << "+" << SYM::get_centering(space_group_string)[i][k];
        vec_pos.push_back(tmp_Wyckoff_equation.str());
        tmp_Wyckoff_equation.str(std::string());
        tmp_Wyckoff_equation.clear();
      }
      equations.push_back(vec_pos);
    }
  }

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << ": Wyckoff positions with centering(s): ";
    for(uint i=0;i<equations.size();i++){
      cerr << aurostd::joinWDelimiter(equations[i],",") << endl;
    }
  }

}
//DX20191029 - END

// ******************************************************************************
// numberOccupiedSitesInConventionalCell() //DX20200601
// ******************************************************************************
namespace SYM {
  uint numberOccupiedSitesInConventionalCell(const vector<wyckoffsite_ITC>& Wyckoff_sites){

    // get the number of occupied sites from the Wyckoff positions
    // i.e., number occupied sites in the conventional cell

    uint natoms = 0;
    for(uint i=0;i<Wyckoff_sites.size();i++){ natoms+=Wyckoff_sites[i].multiplicity; }
    return natoms;

  }
}

// ******************************************************************************
// numberEachTypeFromWyckoff() //DX20200601
// ******************************************************************************
namespace SYM {
  vector<uint> numberEachTypeFromWyckoff(const vector<wyckoffsite_ITC>& Wyckoff_sites){

    // get the number of each atom type from the Wyckoff positions

    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
    compare::groupWyckoffPositions(Wyckoff_sites, grouped_Wyckoff_positions);

    vector<uint> number_of_each_type;
    for(uint i=0;i<grouped_Wyckoff_positions.size();i++){
      uint number_of_type =0;
      for(uint j=0;j<grouped_Wyckoff_positions[i].multiplicities.size();j++){
        number_of_type += grouped_Wyckoff_positions[i].multiplicities[j];
      }
      number_of_each_type.push_back(number_of_type);
    }
    return number_of_each_type;
  }
}

//DX 20191029 - START
// ******************************************************************************
// findWyckoffEquations
// ******************************************************************************
namespace SYM {
  vector<string> findWyckoffEquations(uint space_group_number, string& space_group_setting,
      string& Wyckoff_letter, uint Wyckoff_multiplicity, bool convert2frac){ //DX20200423 - add convert2frac

    SymmetryInformationITC ITC_sym_info;
    ITC_sym_info.initsgs(space_group_setting);
    string spacegroupstring = ITC_sym_info.gl_sgs[space_group_number - 1];
    return findWyckoffEquations(spacegroupstring, Wyckoff_letter, Wyckoff_multiplicity, convert2frac);
  }
}
//DX20191029 END

// ---------------------------------------------------------------------------
namespace SYM {
  vector<string> findWyckoffEquations(string& spacegroupstring, string& Wyckoff_letter, uint Wyckoff_multiplicity, bool convert2frac, bool keep_multiplication_symbol){ //DX20200423 - add convert2frac
    vector<string> positions;
    vector<vector<string> > testvvstring = get_wyckoff_pos(spacegroupstring, Wyckoff_multiplicity, Wyckoff_letter);
    vector<vector<vector<string> > > testvvvstring; testvvvstring.push_back(testvvstring); // quick fix, should probably fix (DX)
    vector<vector<vector<vector<sdouble> > > > testvvvsd = convert_wyckoff_pos_sd(testvvvstring);
    for (uint j = 0; j < testvvvsd.size(); j++) {
      for (uint m = 0; m < testvvvsd[j].size(); m++) {
        string equation = "";
        vector<string> eqn;
        for (uint k = 0; k < 3; k++) { 
          eqn.push_back(SYM::formatWyckoffPosition(testvvvsd[j][m][k], convert2frac, keep_multiplication_symbol)); //DX 20190723 - condensed lines below (in previous version) into function //DX20200423 - add convert2frac, keep_mulitiplication_symbol
        }
        equation = aurostd::joinWDelimiter(eqn,",");
        positions.push_back(equation);
      }
    }
    return positions;
  }
} //namespace SYM
//DX20190128 END

//DX20190723 START
// ******************************************************************************
// formatWyckoffPosition
// ******************************************************************************
namespace SYM {
  string formatWyckoffPosition(const vector<sdouble>& sd_coordinate, bool convert2frac, bool keep_multiplication_symbol,int precision){ //CO20220607 - added precision

    // put variables first (e.g., 1/2+x -> x+1/2)
    // reduce non-variables (e.g., 0.5+0.25 -> 0.75)
    // convert doubles to fractions, if possible (e.g., 0.5 -> 1/2)
    // this function will be circumvented when symbolic math is integrated
    //DX20190723


    stringstream ss_eqn;ss_eqn.precision(precision);  //CO20220607 - added precision
    string coordinate = "";
    vector<string> vec_coord;
    double running_double = 0.0;
    bool double_only = true; //DX20200609 - set default to true and toggle for any variable in j

    for(uint j=0;j<sd_coordinate.size();j++){
      sdouble sd_num = sd_coordinate[j];
      // ---------------------------------------------------------------------------
      // if a variable with positive unit factor (1x -> x) 
      if(sd_num.chr != '\0' && aurostd::abs(sd_num.dbl-1)<_ZERO_TOL_){
        ss_eqn << sd_num.chr;
        vec_coord.push_back(ss_eqn.str());
        double_only = false; //DX20200609
      }
      // ---------------------------------------------------------------------------
      // if a variable with negative unit factor (-1x -> -x)
      else if(sd_num.chr != '\0' && aurostd::abs(sd_num.dbl+1)<_ZERO_TOL_){
        ss_eqn << "-" << sd_num.chr;
        vec_coord.push_back(ss_eqn.str());
        double_only = false; //DX20200609
      }
      // ---------------------------------------------------------------------------
      // if a number (no variable)
      else if(sd_num.chr == '\0'){
        running_double+=sd_num.dbl;
        //DX20200609 [OBSOLETE - doesn't account for other j's] double_only = true;
      }
      // ---------------------------------------------------------------------------
      // if a variable with a scale (2x -> 2x or -0.5x -> -0.5x) 
      else {
        if(keep_multiplication_symbol){ //DX20200428
          ss_eqn << sd_num.dbl << "*" << sd_num.chr;
        }
        else{
          ss_eqn << sd_num.dbl << sd_num.chr;
        }
        vec_coord.push_back(ss_eqn.str());
        double_only = false; //DX20200609
      }
      ss_eqn.str("");
    }

    // ---------------------------------------------------------------------------
    // reduce double, i.e., bring-in-cell; consider positive numbers only, negative numbers are not possible
    // (constants in Wyckoff position are always positive)
    while(running_double>1.0 || aurostd::isequal(running_double,1.0,_ZERO_TOL_)){
      running_double-=1.0;
    }
    //DX20191030 - added case for double only AND not double only - START
    if(double_only){
      if(convert2frac){ //DX20190419
        string running_frac = aurostd::dbl2frac(running_double,false);
        ss_eqn << running_frac;
      } //DX20190419
      //DX20190419 - START
      else{
        ss_eqn << running_double;
      }
      //DX20190419 - END
      //DX20190419 [OBSOLETE] string running_frac = aurostd::dbl2frac(running_double,false); //DX20190724 - now namespace aurostd
      //DX20190419 [OBSOLETE] ss_eqn << running_frac;
      vec_coord.push_back(ss_eqn.str());
      ss_eqn.str("");
    }
    else{
      if(!aurostd::isequal(running_double,0.0,_ZERO_TOL_)){ //DX20190718 - don't add +0 to the end
        string running_frac = aurostd::dbl2frac(running_double,false); //DX20190724 - now namespace aurostd
        ss_eqn << running_frac;
        vec_coord.push_back(ss_eqn.str());
        ss_eqn.str("");
      }
    }
    //DX20191030 - added case for double only AND not double only - END

    coordinate = aurostd::joinWDelimiter(vec_coord,"+");
    // ---------------------------------------------------------------------------
    // clean up cases of -+ or +-
    coordinate = aurostd::StringSubst(coordinate,"-+","-"); 
    coordinate = aurostd::StringSubst(coordinate,"+-","-"); 

    return coordinate;
  }
}

//DX20190708 START
// ******************************************************************************
// reorderWyckoffPosition
// ******************************************************************************
namespace SYM {
  string reorderWyckoffPosition(const string& orig_position){

    // reorder Wyckoff position to a standard format (necessary for string comparisons) 
    // put variables first (e.g., 1/2+x -> x+1/2)
    // reduce non-variables (e.g., 0.5+0.25 -> 0.75)
    // convert doubles to fractions, if possible (e.g., 0.5 -> 1/2)
    // this function will be circumvented when symbolic math is integrated
    //DX20190708

    stringstream message;

    // ---------------------------------------------------------------------------
    // split position via commas 
    vector<string> tokens;
    aurostd::string2tokens(orig_position,tokens,",");

    if(tokens.size()!=3){ 
      message << "Wyckoff position must have 3 fields (e.g., \"x, y, z\"), input: " << orig_position;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
    }

    // ---------------------------------------------------------------------------
    // split position into x,y,z components (via commas)
    string reordered_position = "";
    vector<string> eqn;
    for(uint i=0;i<tokens.size();i++){
      // ---------------------------------------------------------------------------
      // split into equation entities (number, variable) 
      vector<sdouble> sd_coordinate = simplify(tokens[i]);
      eqn.push_back(SYM::formatWyckoffPosition(sd_coordinate)); //DX20190723 - condensed lines below (in previous version) into function
    }
    reordered_position = aurostd::joinWDelimiter(eqn,",");

    return reordered_position;
  }
}
//DX20190708 END

// ******************************************************************************
// shiftWyckoffPositions 
// ******************************************************************************
namespace SYM {
  bool shiftWyckoffPositions(deque<deque<_atom> >& equivalent_atoms_shifted, xvector<double>& previous_shift, 
      xvector<double>& new_shift){
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    // == Subtract shift from previous iterations == //
    for (uint i = 0; i < equivalent_atoms_shifted.size(); i++) {
      for (uint j = 0; j < equivalent_atoms_shifted[i].size(); j++) {
        //DX20190905 [OBSOLETE-no more mod_one_xvec] equivalent_atoms_shifted[i][j].fpos = SYM::mod_one_xvec(equivalent_atoms_shifted[i][j].fpos + previous_shift);
        equivalent_atoms_shifted[i][j].fpos = equivalent_atoms_shifted[i][j].fpos + previous_shift; //DX20190905
        BringInCellInPlace(equivalent_atoms_shifted[i][j].fpos); //DX20190905
      }
    }
    // == Apply new origin shift == //
    for (uint i = 0; i < equivalent_atoms_shifted.size(); i++) {
      for (uint j = 0; j < equivalent_atoms_shifted[i].size(); j++) {
        //DX20190905 [OBSOLETE-no more mod_one_xvec] equivalent_atoms_shifted[i][j].fpos = SYM::mod_one_xvec(equivalent_atoms_shifted[i][j].fpos - new_shift);
        equivalent_atoms_shifted[i][j].fpos = equivalent_atoms_shifted[i][j].fpos - new_shift; //DX20190905
        BringInCellInPlace(equivalent_atoms_shifted[i][j].fpos); //DX20190905
      }
    }
    if(LDEBUG) {
      cerr << "SYM::shiftWyckoffPositions: Shifted atoms: " << endl;
      for(uint i=0;i<equivalent_atoms_shifted.size();i++){
        xb();
        print(equivalent_atoms_shifted[i]);
      }
      cerr << endl;
    } 
    //origin_shift_index++;
    return true;
  }
}

// ******************************************************************************
// findWyckoffPositions 
// ******************************************************************************
// Map equivalent atoms in the xstructure to the ITC Wyckoff positions
namespace SYM {
  bool findWyckoffPositions(xstructure& CCell, deque<_atom>& atomicbasis, vector<vector<vector<string> > >& tmpvvvstring, 
      deque<deque<_atom> >& equivalent_atoms, deque<deque<_atom> >& equivalent_atoms_shifted, 
      bool& foundspacegroup, string& spacegroupstring, bool& orig_origin_shift, xvector<double>& OriginShift, 
      vector<int>& wyckoffmult, vector<string>& wyckoffsymbols, vector<wyckoffsite_ITC>& wyckoffVariables, 
      deque<_atom>& wyckoffPositionsVector, vector<string>& wyckoffSymbols, ostringstream& woss, 
      bool& obverse_force_transformed){

    bool LDEBUG = (FALSE || XHOST.DEBUG);

    bool found_wyckoff = true;
    bool first_wyckoff = true;
    uint found_position = 0;

    // Loop through equivalent atom sets
    for (uint i = 0; i < equivalent_atoms_shifted.size(); i++) {
      int mtmp = equivalent_atoms_shifted[i].size();

      // Find set of Wyckoff positions in the ITC with the same multiplicity 
      vector<int> wyckoffsymbols_mult_index;  //wyckoff symbols with specified multiplicity index
      for (uint ixx = 0; ixx < wyckoffmult.size(); ixx++) {
        if(wyckoffmult[ixx] == mtmp) {
          wyckoffsymbols_mult_index.push_back(ixx);
        }
      }
      first_wyckoff = true;
      uint found_count = 0;
      wyckoffsite_ITC tmp_Wyckoff_site;
      _atom wyckoff_atom;
      found_wyckoff = true;
      tmpvvvstring.clear();

      // Get string of Wyckoff positions
      tmpvvvstring = get_wyckoff_pos(spacegroupstring, mtmp);
      if(tmpvvvstring.size() == 0) {
        if(LDEBUG) { cerr << "SYM::findWyckoffPositions: Could not find wyckoff position based on multiplicity provided [1]. Return string = " << tmpvvvstring.size() << endl; }
        foundspacegroup = false;
        found_wyckoff = false;
        wyckoffVariables.clear();
        wyckoffPositionsVector.clear();
        woss.str("");
        wyckoffSymbols.clear();
        break;
      }
      //xb();
      //print_wyckoff_pos(tmpvvvstring);
      //xb();

      // Convert string to string and doubles to perform math operations on them
      vector<vector<vector<vector<sdouble> > > > tmpvvvsd = convert_wyckoff_pos_sd(tmpvvvstring);
      xvector<double> xyz_params;

      // Loop over the ITC Wyckoff sites with a give multiplicity
      for (uint j = 0; j < tmpvvvsd.size(); j++) {  //LOOP OF WYCKOFF POSITIONS WITH MULTIPLICITY EQUAL TO EQUIVALENT ATOM SET:
        deque<_atom> tmp_equivalent_atoms_shifted = equivalent_atoms_shifted[i];
        vector<vector<vector<sdouble> > > wyckoff_set = tmpvvvsd[j];
        uint variable_changes = 0;
        bool x_variable_set = false;
        bool y_variable_set = false;
        bool z_variable_set = false;

        // Loop over the coordinates for a given ITC Wyckoff set
        for (int m = 0; m < (int)wyckoff_set.size(); m++) {
          bool contains_variable = false;
          bool contains_x_variable = false;
          bool contains_y_variable = false;
          bool contains_z_variable = false;

          // Loop over the equivalent atoms in the structure
          for (uint ix = 0; ix < tmp_equivalent_atoms_shifted.size(); ix++) {  //PERMUTE THROUGH ATOMS IN SET (first may not work)
            xmatrix<double> W(3, 4);
            for (int k = 0; k < 3; k++) {  //TAKE FIRST SYMMETRY-RELATED SITE TO ASSIGN VARIABLES (HENCE: [j][0])
              W(k + 1, 4) = tmp_equivalent_atoms_shifted[ix].fpos(k + 1);
              for (uint l = 0; l < wyckoff_set[m][k].size(); l++) {
                //cerr <<"is " << wyckoff_set[m][k][l].dbl << " " << wyckoff_set[m][k][l].chr << " equal to: " << tmp_equivalent_atoms_shifted[ix].fpos(k+1) << endl;//DEBUG MODE
                // === Compare constent components in Wyckoff positions first === //
                if(wyckoff_set[m][k][l].chr == '\0') {
                  W(k + 1, 4) -= wyckoff_set[m][k][l].dbl;
                } 
                else {
                  int neg_pos = 1;
                  //If the coefficient is negative then multiply both sides by negative one and take RHS modulo 1
                  //cerr << "wyckoff_set[m][k][l].chr: " << wyckoff_set[m][k][l].chr << endl;
                  if(wyckoff_set[m][k][l].chr == 'x') {
                    contains_variable = true;
                    contains_x_variable = true;
                    W(k + 1, 1) = neg_pos * wyckoff_set[m][k][l].dbl;
                  }
                  if(wyckoff_set[m][k][l].chr == 'y') {
                    contains_variable = true;
                    contains_y_variable = true;
                    W(k + 1, 2) = neg_pos * wyckoff_set[m][k][l].dbl;
                  }
                  if(wyckoff_set[m][k][l].chr == 'z') {
                    contains_variable = true;
                    contains_z_variable = true;
                    W(k + 1, 3) = neg_pos * wyckoff_set[m][k][l].dbl;
                  }
                  //DX20190905 [OBSOLETE-no more mod_one_xvec] W(k + 1, 4) = SYM::mod_one(neg_pos * tmp_equivalent_atoms_shifted[ix].fpos(k + 1));
                  W(k + 1, 4) = neg_pos * tmp_equivalent_atoms_shifted[ix].fpos(k + 1); //DX20190905
                  BringInCellInPlace(W(k + 1, 4)); //DX20190905
                }
              }
            }
            // ===== Linear algebra problem: Solve for Wyckoff position (if positions are not constant) ===== //
            vector<xvector<double> > LHS;
            xvector<double> tmp1; tmp1(1) = W(1, 1); tmp1(2) = W(1, 2); tmp1(3) = W(1, 3);
            xvector<double> tmp2; tmp2(1) = W(2, 1); tmp2(2) = W(2, 2); tmp2(3) = W(2, 3);
            xvector<double> tmp3; tmp3(1) = W(3, 1); tmp3(2) = W(3, 2); tmp3(3) = W(3, 3);
            LHS.push_back(tmp1);
            LHS.push_back(tmp2);
            LHS.push_back(tmp3);
            for (uint i = 0; i < 3; i++) {
              for (uint j = 0; j < 3; j++) {
                for (uint k = 0; k < 3; k++) {
                  vector<double> dS_mod;
                  dS_mod.push_back(W(1, 4));
                  dS_mod.push_back(W(2, 4));
                  dS_mod.push_back(W(3, 4));
                  dS_mod[0] = dS_mod[0] + double(i);
                  dS_mod[1] = dS_mod[1] + double(j);
                  dS_mod[2] = dS_mod[2] + double(k);
                  xmatrix<double> W_tmp = xvec2xmat(LHS, dS_mod);
                  //cerr <<"W_tmp: " <<  endl; //DEBUG MODE
                  //cerr << W_tmp << endl; //DEBUG MODE
                  //xb(); //DEBUG MODE
                  ReducedRowEchelonForm(W_tmp,CCell.sym_eps); //DX20190215
                  //cerr << "RR FORM: " << endl; //DEBUG MODE
                  //cerr << W_tmp << endl; //DEBUG MODE
                  xvector<double> solution; solution(1) = (W_tmp(1, 4)); solution(2) = (W_tmp(2, 4)); solution(3) = (W_tmp(3, 4));
                  solution = SYM::minimizeDistanceFractionalMethod(solution); //DX20190613
                  //DX20190613 [OBSOLETE] SYM::PBC(solution); 
                  xvector<double> cart_solution = trasp(CCell.lattice)*solution;

                  xvector<double> ideal_wyckoff = solution;
                  found_wyckoff = true;
                  if((aurostd::abs(W_tmp(1, 1)) < _ZERO_TOL_ && aurostd::abs(W_tmp(1, 2)) < _ZERO_TOL_ && aurostd::abs(W_tmp(1, 3)) < _ZERO_TOL_)){
                    ideal_wyckoff(1) = 0.0;
                  }
                  if((aurostd::abs(W_tmp(2, 1)) < _ZERO_TOL_ && aurostd::abs(W_tmp(2, 2)) < _ZERO_TOL_ && aurostd::abs(W_tmp(2, 3)) < _ZERO_TOL_)){
                    ideal_wyckoff(2) = 0.0;
                  }
                  if((aurostd::abs(W_tmp(3, 1)) < _ZERO_TOL_ && aurostd::abs(W_tmp(3, 2)) < _ZERO_TOL_ && aurostd::abs(W_tmp(3, 3)) < _ZERO_TOL_)){
                    ideal_wyckoff(3) = 0.0;
                  }
                  ideal_wyckoff = SYM::minimizeDistanceFractionalMethod(ideal_wyckoff); //DX20190613
                  //DX20190613 [OBSOLETE] SYM::PBC(ideal_wyckoff);
                  xvector<double> ideal_cart_position = trasp(CCell.lattice)*ideal_wyckoff;
                  if(aurostd::modulus(cart_solution-ideal_cart_position)>CCell.sym_eps){ //DX20190215
                    found_wyckoff = false;
                  }
                  if(found_wyckoff) {
                    //cerr << "FOUND: " << W_tmp << endl;
                    W = W_tmp;
                    break;
                  }
                }
                if(found_wyckoff) {
                  break;
                }
              }
              if(found_wyckoff) {
                break;
              }
            }
            //cerr << "FOUND WYCKOFF: " << found_wyckoff << endl;
            if(found_wyckoff == true) {
              for (int ixx = 1; ixx < 4; ixx++) {
                for (int jx = 1; jx < 4; jx++) {
                  if(aurostd::abs(W(ixx, jx) - 1.0) < _ZERO_TOL_) {
                    if(jx == 1) {
                      if(contains_x_variable && !x_variable_set){
                        xyz_params(1) = W(ixx,4);
                        x_variable_set = true;
                      }
                      if(first_wyckoff == true){
                        tmp_Wyckoff_site.coord(1) = W(ixx, 4);
                      }
                    }
                    if(jx == 2) {
                      //cerr << "y = " << W(ix,4) << " ";
                      if(contains_y_variable && !y_variable_set){
                        xyz_params(2) = W(ixx,4);
                        y_variable_set = true;
                      }
                      if(first_wyckoff == true){
                        tmp_Wyckoff_site.coord(2) = W(ixx, 4);
                      }
                    }
                    if(jx == 3) {
                      if(contains_z_variable && !z_variable_set){
                        xyz_params(3) = W(ixx,4);
                        z_variable_set = true;
                      }
                      if(first_wyckoff == true){
                        tmp_Wyckoff_site.coord(3) = W(ixx, 4);
                      }
                    }
                  }
                }
                if(first_wyckoff == true){
                  //cerr << "FOUND LINEAR SOLUTION: " << endl;
                  //cerr << endl << W << endl;
                  //cerr << endl;
                  first_wyckoff = false;
                  // ========== Store Wyckoff Position ========== //
                  tmp_Wyckoff_site.coord = tmp_equivalent_atoms_shifted[ix].fpos; //DX20171212 - need to updated coord to include non-parametrized Wyckoff positions
                  tmp_Wyckoff_site.index = tmp_equivalent_atoms_shifted[ix].type; //DX20200427
                  tmp_Wyckoff_site.type = tmp_equivalent_atoms_shifted[ix].name;
                  //DX20191010 - update partial occupation value - START
                  if(CCell.partial_occupation_flag){
                    tmp_Wyckoff_site.site_occupation=tmp_equivalent_atoms_shifted[ix].partial_occupation_value;
                  }
                  else{ tmp_Wyckoff_site.site_occupation=1.0; }
                  //DX20191010 - update partial occupation value - END
                  tmp_Wyckoff_site.wyckoffSymbol = wyckoffsymbols[wyckoffsymbols_mult_index[j] - 1];
                  tmp_Wyckoff_site.equations.clear(); //DX20190130 - clear for safety
                  //DX20190128 - add multiplicity, letter, and site symmetry - START 
                  vector<string> tokens,pos_tokens;
                  if(aurostd::string2tokens(tmp_Wyckoff_site.wyckoffSymbol,tokens," ") == 3){
                    tmp_Wyckoff_site.multiplicity = aurostd::string2utype<int>(tokens[0]);
                    tmp_Wyckoff_site.letter = tokens[1];
                    tmp_Wyckoff_site.site_symmetry = tokens[2];
                    vector<string> positions = SYM::findWyckoffEquations(spacegroupstring, tmp_Wyckoff_site.letter, 
                        tmp_Wyckoff_site.multiplicity);
                    for(uint p=0;p<positions.size();p++){
                      aurostd::string2tokens(positions[p],pos_tokens,",");
                      tmp_Wyckoff_site.equations.push_back(pos_tokens);
                    }
                  }
                  //DX20190128 - add multiplicity, letter, and site symmetry - END
                  wyckoff_atom = tmp_equivalent_atoms_shifted[ix];
                }
              }
              if(contains_variable){
                xvector<double> diff;
                //cerr << "contains variable" << endl;
                for (int ixx = 1; ixx < 4; ixx++) {
                  for (int jx = 1; jx < 4; jx++) {
                    if(aurostd::abs(W(ixx, jx) - 1.0) < CCell.sym_eps) { //DX20190215
                      if(contains_x_variable && x_variable_set && jx ==1){
                        diff(1) = xyz_params(1)-W(ixx,4);
                      }
                      if(contains_y_variable && y_variable_set && jx == 2){
                        diff(2) = xyz_params(2)-W(ixx,4);
                      }
                      if(contains_z_variable && z_variable_set && jx == 3){
                        diff(3) = xyz_params(3)-W(ixx,4);
                      }
                    }
                    if(!found_wyckoff){
                      break;
                    }
                  }
                }
                diff = SYM::minimizeDistanceFractionalMethod(diff); //DX20190613
                //DX20190613 [OBSOLETE] SYM::PBC(diff);
                if(aurostd::modulus(trasp(CCell.lattice)*diff)>CCell.sym_eps){ //DX20190215
                  found_wyckoff = false;
                }
              }
              if(!found_wyckoff){
                continue;
              }
              found_count++;
              tmp_equivalent_atoms_shifted.erase(tmp_equivalent_atoms_shifted.begin() + ix);
              wyckoff_set.erase(wyckoff_set.begin() + m);
              m--;
              break;
            }
          }  //for ix: loop through equivalent atoms
          if(found_wyckoff == false) {
            //cerr << "COULDN'T FIND MATCH: equivalent atoms : " << i << " with " << m << endl;
            first_wyckoff = true;
            found_count = 0;
            //found_position = 0;
            // Perhaps we found the wrong parameter for a variable, rearrange equivalent atom order to find new parameter
            if(contains_variable && variable_changes < equivalent_atoms_shifted[i].size() - 1) {
              variable_changes += 1;
              x_variable_set = false;
              y_variable_set = false;
              z_variable_set = false;
              contains_x_variable = false;
              contains_y_variable = false;
              contains_z_variable = false;
              xyz_params.clear();
              tmp_equivalent_atoms_shifted = equivalent_atoms_shifted[i];
              wyckoff_set = tmpvvvsd[j];
              for (uint s = 0; s < variable_changes; s++) {
                tmp_equivalent_atoms_shifted.push_back(tmp_equivalent_atoms_shifted[0]);
                tmp_equivalent_atoms_shifted.erase(tmp_equivalent_atoms_shifted.begin() + 0);
              }
              m = -1;
            }
            else {
              break;
            }
          }
        }  // for m: loop through positions in a given Wyckoff set
        if(found_count == equivalent_atoms_shifted[i].size()) {
          wyckoffVariables.push_back(tmp_Wyckoff_site);
          if(LDEBUG) {
            cerr << "SYM::findWyckoffPosition: Found a Wyckoff position: " << endl << tmp_Wyckoff_site << endl;
          } 
          //cerr << "wyckoff site: ";
          //woss << equivalent_atoms_shifted[i][ix].coord <<"   Atom" <<  equivalent_atoms_shifted[i][ix].type << wyckoffsymbols[wyckoffsymbols_mult_index[j]-1] << endl; //GET WYCKOFF SITE SYMBOL FOR APPROPRIATE MULT

          //woss << equivalent_atoms_shifted[i][ix].coord <<" " <<  CCell.chemical_labels[atoi(equivalent_atoms_shifted[i][ix].name.c_str())] << wyckoffsymbols[wyckoffsymbols_mult_index[j]-1] << endl; //GET WYCKOFF SITE SYMBOL FOR APPROPRIATE MULT
          wyckoffPositionsVector.push_back(wyckoff_atom);                                                                                                                                                   //Modified (RHT)
          woss << setprecision(14) << fixed << "   " << wyckoff_atom.fpos(1) << "   " << wyckoff_atom.fpos(2) << "   " << wyckoff_atom.fpos(3) << "   " << wyckoff_atom.name << tmp_Wyckoff_site.wyckoffSymbol << endl;  //GET WYCKOFF SITE SYMBOL FOR APPROPRIATE MULT
          wyckoffSymbols.push_back(tmp_Wyckoff_site.wyckoffSymbol);
          found_position++;
        }
        if(found_count == equivalent_atoms_shifted[i].size()) {
          break;
        }
      }  //for j : loop through Wyckoff sets with multiplicity
      if(found_count != equivalent_atoms_shifted[i].size() && orig_origin_shift == true) {
        //if(LDEBUG) { cerr << "SYM::findWyckoffPositions WARNING: Equivalent atoms were not matched with ITC Wyckoff positions." << endl; }
        if(LDEBUG) {
          cerr << "SYM::findWyckoffPositions WARNING: Could not match the " << equivalent_atoms_shifted[i][0].name << " atom(s) with multiplicity " << equivalent_atoms_shifted[i].size() << "." << endl;
        }
        first_wyckoff = true;
        wyckoffVariables.clear();
        wyckoffPositionsVector.clear();
        woss.str("");
        wyckoffSymbols.clear();
        if(CCell.lattice_label_ITC == 'R' && obverse_force_transformed == false) {
          transformToObverse(CCell.lattice, equivalent_atoms_shifted);
          obverse_force_transformed = true;
          xvector<double> ReverseOriginShift = -OriginShift;
          equivalent_atoms = shiftSymmetryEquivalentAtoms(equivalent_atoms_shifted, CCell.lattice, ReverseOriginShift, CCell.dist_nn_min, CCell.sym_eps); //DX20190215 - added CCell.sym_eps
          deque<_atom> atoms_tmp;
          for (uint k = 0; k < equivalent_atoms.size(); k++) {
            for (uint l = 0; l < equivalent_atoms[k].size(); l++) {
              atoms_tmp.push_back(equivalent_atoms[k][l]);
            }
          }
          atomicbasis = atoms_tmp;
          CCell.atoms = atoms_tmp;
          i = -1;
          found_position=0; //DX20170912
        } 
        else {
          break;
        }
      } 
      else if(found_count != equivalent_atoms_shifted[i].size() && orig_origin_shift == false) {
        first_wyckoff = true;
        wyckoffVariables.clear();
        wyckoffPositionsVector.clear();
        woss.str("");
        wyckoffSymbols.clear();
        //sum_wyckoff_letters.push_back(100);
        break;
      }
      if(tmpvvvstring.size() == 0) {
        if(LDEBUG) { cerr << "SYM::findWyckoffPositions: Could not find wyckoff position based on multiplicity provided [2]. Return string = " << tmpvvvstring.size() << endl; }
        break;
      }
    }
    return (found_position == equivalent_atoms_shifted.size());
  }
} //namespace SYM

// ******************************************************************************
// GetSameSymmetryWyckoffLetters
// ******************************************************************************
namespace SYM {
  vector<vector<vector<string> > > GetSameSymmetryWyckoffLetters(uint space_group_number, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, uint setting){;
    // get possible Wyckoff letters
    stringstream axis_cell;
    axis_cell.str(std::string());
    //DX20180806 [OBSOLETE] axis_cell << CCell.lattice_label_ITC << cell_choice;
    axis_cell << setting; //DX20180806 - use setting
    SymmetryInformationITC ITC_sym_info; //DX20190215
    ITC_sym_info.initsgs(axis_cell.str()); //DX20190215

    //DX20190215 [OBSOLETE] using SYM::gl_sgs;
    string spacegroupstring = ITC_sym_info.gl_sgs[space_group_number - 1];
    vector<vector<vector<string> > > grouped_possible_Wyckoff_letters;
    for(uint i=0;i<grouped_Wyckoff_positions.size();i++){
      vector<vector<string> > possible_Wyckoff_letters;
      vector<vector<string> > grouped_possible_Wyckoff_site_symmetries;
      for(uint j=0;j<grouped_Wyckoff_positions[i].site_symmetries.size();j++){
        vector<string> site_syms;
        vector<string> letters;
        vector<string> positions;
        SYM::get_certain_wyckoff_pos(spacegroupstring, grouped_Wyckoff_positions[i].multiplicities[j], grouped_Wyckoff_positions[i].site_symmetries[j], site_syms, letters, positions);
        possible_Wyckoff_letters.push_back(letters);
      }
      grouped_possible_Wyckoff_letters.push_back(possible_Wyckoff_letters);
    }
    for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){
      for(uint j=0;j<grouped_possible_Wyckoff_letters[i].size();j++){
        ::print(grouped_possible_Wyckoff_letters[i][j]);
      }
    }
    return grouped_possible_Wyckoff_letters; 
  }
}

//// ******************************************************************************
//// GetSameSymmetryWyckoffLetters
//// ******************************************************************************
//namespace SYM {
//vector<vector<vector<string> > > GetSameSymmetryWyckoffLetters(uint space_group_number, string Wyckoff_letter_string, uint setting){;
//vector<vector<string> > all_Wyckoff_sets;
//vector<string> tokens, per_atom_tokens;
//aurostd::string2tokens(Wyckoff_letter_string, tokens, ";");
//
//for(uint i=0;i<tokens.size();i++){
//vector<string> Wyckoff_set;
//aurostd::string2tokens(tokens[i],per_atom_tokens,",");
//for(uint j=0;j<per_atom_tokens.size();j++){
//Wyckoff_set.push_back(per_atom_tokens[j]);
//}
//all_Wyckoff_sets.push_back(Wyckoff_set);
//}
//
//// get possible Wyckoff letters
//stringstream axis_cell;
//axis_cell.str(std::string());
////DX20180806 [OBSOLETE] axis_cell << CCell.lattice_label_ITC << cell_choice;
//axis_cell << setting; //DX20180806 - use setting
//SYM::initsgs(axis_cell.str());
////extern vector<string> gl_sgs;
//using SYM::gl_sgs;
//string spacegroupstring = gl_sgs[space_group_number - 1];
//vector<vector<vector<string> > > grouped_possible_Wyckoff_letters;
//for(uint i=0;i<grouped_Wyckoff_positions.size();i++){
//vector<vector<string> > possible_Wyckoff_letters;
//vector<vector<string> > grouped_possible_Wyckoff_site_symmetries;
//for(uint j=0;j<grouped_Wyckoff_positions[i].site_symmetries.size();j++){
//vector<string> site_syms;
//vector<string> letters;
//vector<string> positions;
//SYM::get_certain_wyckoff_pos(spacegroupstring, grouped_Wyckoff_positions[i].multiplicities[j], grouped_Wyckoff_positions[i].site_symmetries[j], site_syms, letters, positions);
//possible_Wyckoff_letters.push_back(letters);
//}
//grouped_possible_Wyckoff_letters.push_back(possible_Wyckoff_letters);
//}
//for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){
//for(uint j=0;j<grouped_possible_Wyckoff_letters[i].size();j++){
//::print(grouped_possible_Wyckoff_letters[i][j]);
//}
//}
//return grouped_possible_Wyckoff_letters; 
//}
//}

// ******************************************************************************
// blank
// ******************************************************************************
namespace SYM {
  bool blank(string str_in) {
    bool out = true;
    for (uint i = 0; i < str_in.size(); i++) {
      if(str_in[i] != ' ')
        out = false;
    }
    return out;
  }
} //namespace SYM

// ******************************************************************************
// containschar
// ******************************************************************************
//If string contains alphabetic character
namespace SYM {
  bool containschar(string str_in) {
    bool contains = false;
    for (uint i = 0; i < str_in.length(); i++) {
      if(isalpha(str_in[i]))
        contains = true;
    }
    return contains;
  }
} //namespace SYM

//DX20210421 [OBSOLETE]// ******************************************************************************
//DX20210421 [OBSOLETE]// intinvec
//DX20210421 [OBSOLETE]// ******************************************************************************
//DX20210421 [OBSOLETE]namespace SYM {
//DX20210421 [OBSOLETE]  bool intinvec(vector<int> vint, int n) {
//DX20210421 [OBSOLETE]    bool contains = false;
//DX20210421 [OBSOLETE]    for (uint i = 0; i < vint.size(); i++) {
//DX20210421 [OBSOLETE]      if(vint[i] == n)
//DX20210421 [OBSOLETE]        contains = true;
//DX20210421 [OBSOLETE]    }
//DX20210421 [OBSOLETE]    return contains;
//DX20210421 [OBSOLETE]  }
//DX20210421 [OBSOLETE]} //namespace SYM

// ******************************************************************************
// whichchar
// ******************************************************************************
namespace SYM {
  char whichchar(const string& str_in) {
    char c = '\0';
    for (uint i = 0; i < str_in.length(); i++) {
      if(isalpha(str_in[i]))
        c = str_in[i];
    }
    return c;
  }
} //namespace SYM

// ******************************************************************************
// isop
// ******************************************************************************
namespace SYM {
  bool isop(char c) {
    bool isanoperator = false;

    //char plus='+';

    if(c == '+' || c == '-' || c == '*' || c == '/')
      isanoperator = true;

    return isanoperator;
  }
} //namespace SYM

// ******************************************************************************
// whichnum
// ******************************************************************************
namespace SYM {
  double whichnum(string str_in) {  //returns the double that is contained in string. e.g., .25x returns .25
    bool neg = false;
    if(str_in[0] == '-') {
      neg = true;
    }
    double num = 0;
    ostringstream oss;
    int numcount = 0;
    for (uint i = 0; i < str_in.length(); i++) {
      if(!isalpha(str_in[i]) && !isop(str_in[i])) {
        oss << str_in[i];
        numcount++;
      }
    }
    num = atof(oss.str().c_str());
    if(numcount == 0)
      num = 1;
    if(neg == true)
      num = -num;
    return num;
  }
} //namespace SYM

// ******************************************************************************
// havechar
// ******************************************************************************
//if the character c is contained in the string
namespace SYM {
  bool havechar(string str_in, char c) {
    bool contains = false;
    for (uint i = 0; i < str_in.length(); i++) {
      if(str_in[i] == c)
        contains = true;
    }
    return contains;
  }
} //namespace SYM

// ******************************************************************************
// convert_wyckoff_pos
// ******************************************************************************
//For now is void, but will be matrix
namespace SYM {
  void convert_wyckoff_pos(vector<vector<vector<string> > > wyckoff_positions) {
    cerr << "*********************************" << endl;
    ostringstream oss;
    //  atof(tmpstr.c_str())
    for (uint k = 0; k < wyckoff_positions.size(); k++) {
      for (uint j = 0; j < wyckoff_positions[k].size(); j++) {
        for (uint i = 0; i < wyckoff_positions[k][j].size(); i++) {
          oss << wyckoff_positions[k][j][i];
        }
        string temp = oss.str();
        if(!containschar(temp))
          cerr << atof(temp.c_str());
        else {
          vector<string> tempvec = splitstring(temp);
          for (uint i = 0; i < tempvec.size(); i++) {
            if(!isalpha(tempvec[i][0]))
              cerr << atof(tempvec[i].c_str());
            else
              cerr << tempvec[i];
          }
        }
        oss.str("");
        cerr << " ";
      }
      cerr << endl;
    }
  }
} //namespace SYM

// ******************************************************************************
// splitstring
// ******************************************************************************
namespace SYM {
  vector<string> splitstring(string str) {
    ostringstream oss;
    vector<string> temp;
    for (uint i = 0; i < str.length(); i++) {
      oss << str[i];
      //if(str[i+1]=='-' || str[i+1]=='+' || str[i+1]=='\0' || str[i+1]==' ')
      if(!blank(oss.str()))
      { //CO20200106 - patching for auto-indenting
        //DX20200924 [OBSOLETE] if(str[i + 1] == '\0' || str[i + 1] == ' ' || str[i + 1] == '+' || str[i + 1] == '-')
        //DX20200924 - START
        // distinguish between empty lines and spaces
        if(str[i + 1] == '\0' || str[i + 1] == ' '){
          temp.push_back(oss.str());
          oss.str("");
        }
        // check for + or -
        else if(str[i + 1] == '+' || str[i + 1] == '-') {
          // ensure their isn't a multiplication/division prior to the + or -
          // (otherwise this indicates the sign of the value)
          if(str[i] != '*' && str[i] != '/'){
            temp.push_back(oss.str());
            oss.str("");
          }
        }
        //DX20200924 - END
      }
    }

    //if(str.length()==1 && str[0]=='0')
    //temp.push_back("0");

    return temp;
  }
} //namespace SYM

// ******************************************************************************
// splitstring_2
// ******************************************************************************
namespace SYM {
  vector<string> splitstring_2(string str) {
    ostringstream oss;
    vector<string> temp;
    for (uint i = 0; i < str.length(); i++) {
      oss << str[i];
      //if(str[i+1]=='-' || str[i+1]=='+' || str[i+1]=='\0' || str[i+1]==' ')
      if(str[i + 1] == '\0' || str[i + 1] == ' ')
      { //CO20200106 - patching for auto-indenting
        temp.push_back(oss.str());
        //cerr << oss.str();// << endl;
        oss.str("");
      }
    }

    //if(str.length()==1 && str[0]=='0')
    //temp.push_back("0");

    return temp;
  }
} //namespace SYM

//DX20190724 [MOVED TO AUROSTD] - dbl2frac() Double to Fraction (Overloaded)

//DX20200313 [reformatted and moved to AUROSTD]// ******************************************************************************
//DX20200313 [reformatted and moved to AUROSTD]// dbl2frac Double to Fraction (Overloaded)
//DX20200313 [reformatted and moved to AUROSTD]// ******************************************************************************
//DX20200313 [reformatted and moved to AUROSTD]//NEED TO REWRITE/check FRAC2DBL
//DX20200313 [reformatted and moved to AUROSTD]namespace SYM {
//DX20200313 [reformatted and moved to AUROSTD]  double frac2dbl(string str) {  //string in the form "a/b" with minus if negative
//DX20200313 [reformatted and moved to AUROSTD]    //  if(str[0]==' '){
//DX20200313 [reformatted and moved to AUROSTD]    //str.erase(str.begin(),str.begin()+1);
//DX20200313 [reformatted and moved to AUROSTD]    //}
//DX20200313 [reformatted and moved to AUROSTD]
//DX20200313 [reformatted and moved to AUROSTD]    double out = 0;
//DX20200313 [reformatted and moved to AUROSTD]    bool neg = false;
//DX20200313 [reformatted and moved to AUROSTD]    stringstream ss_num, ss_den;
//DX20200313 [reformatted and moved to AUROSTD]    cleanupstring(str);
//DX20200313 [reformatted and moved to AUROSTD]    if(str[0] == '-') {
//DX20200313 [reformatted and moved to AUROSTD]      neg = true;
//DX20200313 [reformatted and moved to AUROSTD]      str.erase(str.begin(), str.begin() + 1);
//DX20200313 [reformatted and moved to AUROSTD]    }
//DX20200313 [reformatted and moved to AUROSTD]    if(str[0] == '+') {
//DX20200313 [reformatted and moved to AUROSTD]      str.erase(str.begin(), str.begin() + 1);
//DX20200313 [reformatted and moved to AUROSTD]    }
//DX20200313 [reformatted and moved to AUROSTD]    if(havechar(str, '/')) {
//DX20200313 [reformatted and moved to AUROSTD]      double numerator = 0;
//DX20200313 [reformatted and moved to AUROSTD]      double denominator = 0;
//DX20200313 [reformatted and moved to AUROSTD]      uint slash = whereischar(str, '/');
//DX20200313 [reformatted and moved to AUROSTD]      for (uint i = 0; i < slash; i++)
//DX20200313 [reformatted and moved to AUROSTD]        ss_num << str[i];
//DX20200313 [reformatted and moved to AUROSTD]      for (uint i = 0; i < (str.size() - slash); i++)
//DX20200313 [reformatted and moved to AUROSTD]        ss_den << str[i + slash + 1];
//DX20200313 [reformatted and moved to AUROSTD]      string num_string = ss_num.str();
//DX20200313 [reformatted and moved to AUROSTD]      const char* num = num_string.c_str();
//DX20200313 [reformatted and moved to AUROSTD]      numerator = atof(num);
//DX20200313 [reformatted and moved to AUROSTD]      string den_string = ss_den.str();
//DX20200313 [reformatted and moved to AUROSTD]      const char* den = den_string.c_str();
//DX20200313 [reformatted and moved to AUROSTD]      denominator = atof(den);
//DX20200313 [reformatted and moved to AUROSTD]      if(neg == true) {
//DX20200313 [reformatted and moved to AUROSTD]        out = -numerator / denominator;
//DX20200313 [reformatted and moved to AUROSTD]      } else {
//DX20200313 [reformatted and moved to AUROSTD]        out = numerator / denominator;
//DX20200313 [reformatted and moved to AUROSTD]      }
//DX20200313 [reformatted and moved to AUROSTD]    } else {
//DX20200313 [reformatted and moved to AUROSTD]      char dbl[256];
//DX20200313 [reformatted and moved to AUROSTD]      for (uint i = 0; i < str.size(); i++)
//DX20200313 [reformatted and moved to AUROSTD]        dbl[i] = str[i];
//DX20200313 [reformatted and moved to AUROSTD]      out = atof(dbl);
//DX20200313 [reformatted and moved to AUROSTD]    }
//DX20200313 [reformatted and moved to AUROSTD]    return out;
//DX20200313 [reformatted and moved to AUROSTD]  }
//DX20200313 [reformatted and moved to AUROSTD]}

//DX20210422 [OBSOLETE]// ******************************************************************************
//DX20210422 [OBSOLETE]// distance_between_points
//DX20210422 [OBSOLETE]// ******************************************************************************
//DX20210422 [OBSOLETE]namespace SYM {
//DX20210422 [OBSOLETE]  double distance_between_points(const xvector<double>& a, const xvector<double>& b) {
//DX20210422 [OBSOLETE]    double dist = 0;
//DX20210422 [OBSOLETE]    if(a.urows != b.urows) {
//DX20210422 [OBSOLETE]      throw aurostd::xerror(__AFLOW_FILE__,XPID+"SYM::distance_between_points():","Vectors are not equal length.",_VALUE_RANGE_);
//DX20210422 [OBSOLETE]    }
//DX20210422 [OBSOLETE]    for (uint i = 1; i <= (uint)a.urows; i++) {
//DX20210422 [OBSOLETE]      dist += (a(i) - b(i)) * (a(i) - b(i));
//DX20210422 [OBSOLETE]    }
//DX20210422 [OBSOLETE]    dist = sqrt(dist);
//DX20210422 [OBSOLETE]    return dist;
//DX20210422 [OBSOLETE]  }
//DX20210422 [OBSOLETE]} //namespace SYM

//DX20210422 [OBSOLETE]// ******************************************************************************
//DX20210422 [OBSOLETE]// xstring (Overloaded)
//DX20210422 [OBSOLETE]// ******************************************************************************
//DX20210422 [OBSOLETE]namespace SYM {
//DX20210422 [OBSOLETE]  void xstring(ostream& output, xmatrix<double> a) {
//DX20210422 [OBSOLETE]    for (uint i = 1; i <= (uint)a.urows; i++) {
//DX20210422 [OBSOLETE]      for (uint j = 1; j <= (uint)a.lrows; j++) {
//DX20210422 [OBSOLETE]        if(aurostd::abs(a(i, j)) < 1e-8) { a(i, j) = 0; }
//DX20210422 [OBSOLETE]        output << std::showpoint << setprecision(14) << a(i, j) << " ";
//DX20210422 [OBSOLETE]      }
//DX20210422 [OBSOLETE]      cerr << endl;
//DX20210422 [OBSOLETE]    }
//DX20210422 [OBSOLETE]  }
//DX20210422 [OBSOLETE]} //namespace SYM

// ******************************************************************************
// minus_one
// ******************************************************************************
namespace SYM {
  void minus_one(deque<_atom>& atomicBasis, int lvec) {
    //DX double tol = 1e-6;
    double tol = 1e-10;
    for (uint i = 0; i < atomicBasis.size(); i++) {
      if(atomicBasis[i].fpos(lvec) > tol) {
        atomicBasis[i].fpos(lvec) = 1 - atomicBasis[i].fpos(lvec);
      }
    }
  }
} //namespace SYM

// ******************************************************************************
// swap_columns
// ******************************************************************************
namespace SYM {
  void swap_columns(deque<_atom>& atomicBasis, int col1, int col2) {
    double tmp = 0;
    for (uint i = 0; i < atomicBasis.size(); i++) {
      tmp = atomicBasis[i].fpos(col1);
      atomicBasis[i].fpos(col1) = atomicBasis[i].fpos(col2);
      atomicBasis[i].fpos(col2) = tmp;
    }
  }
} //namespace SYM

// ******************************************************************************
// rearrange_columns
// ******************************************************************************
namespace SYM {
  void rearrange_columns(deque<_atom>& atomicBasis, int c1, int c2, int c3) {
    double tmp1 = 0;
    double tmp2 = 0;
    double tmp3 = 0;
    for (uint i = 0; i < atomicBasis.size(); i++) {
      tmp1 = atomicBasis[i].fpos(c1 + 1);
      tmp2 = atomicBasis[i].fpos(c2 + 1);
      tmp3 = atomicBasis[i].fpos(c3 + 1);
      atomicBasis[i].fpos(1) = tmp1;
      atomicBasis[i].fpos(2) = tmp2;
      atomicBasis[i].fpos(3) = tmp3;
    }
  }
} //namespace SYM

// ******************************************************************************
// simplify
// ******************************************************************************
namespace SYM{
  vector<sdouble> simplify(const string& _str) {
    string str=_str;
    if(str[0] == ' ') {
      str.erase(str.begin(), str.begin() + 1);
    }
    // if contains double negative, make positive //DX20200622
    if(aurostd::substring2bool(str,"--")){
      aurostd::StringSubst(str,"--","+");
    }
    sdouble out;
    vector<sdouble> out_vec;
    double num = 0;
    vector<string> splitstr = splitstring(str);
    for (uint S = 0; S < splitstr.size(); S++) {
      num = 0;  //creates some problem for unreduced strings (e.g, .5+.3)
      out.chr = '\0';
      bool skip = false;
      ostringstream oss;
      for (uint i = 0; i < splitstr[S].length(); i++) {
        oss << splitstr[S][i];
      }
      if(!containschar(oss.str())) {
        if(havechar(oss.str(), '/')) {
          //DX20200313 [OBSOLETE] num += frac2dbl(oss.str());
          num += aurostd::frac2dbl(oss.str()); // DX20200313
        } else if(havechar(oss.str(), '*')) {
          vector<string> tokens; aurostd::string2tokens(oss.str(),tokens,"*");
          if(tokens.size()==2){ num += atof(tokens[0].c_str())*atof(tokens[1].c_str()); }
          else { throw aurostd::xerror(__AFLOW_FILE__,"SYM::simplify","Expected two fields surrounding '*': "+oss.str(),_RUNTIME_ERROR_); }
        } else {
          num += atof(oss.str().c_str());
        }
        oss.str("");
      } else {
        skip = true;
        if(havechar(oss.str(), '/')) {
          //DX20200313 [OBSOLETE] out.dbl = frac2dbl(oss.str());
          out.dbl = aurostd::frac2dbl(oss.str()); //DX20200313
        } else {
          out.dbl = whichnum(oss.str());
        }
        out.chr = whichchar(oss.str());
        oss.str("");
      }
      if(skip == false) {
        out.dbl = num;
      }
      if(aurostd::abs(out.dbl) > .00001) {
        out_vec.push_back(out);
      }
    }
    if(out_vec.size() == 0) {
      out.dbl = 0;
      out_vec.push_back(out);
    }
    return out_vec;
  }
} //namespace SYM

// ******************************************************************************
// print_wyckoff_pos
// ******************************************************************************
namespace SYM {
  void print_wyckoff_pos(vector<vector<vector<string> > > wyckoff_positions) {
    //Wyckoff_positions[tuplets][singlets][singlet characters]
    for (uint k = 0; k < wyckoff_positions.size(); k++) {
      for (uint j = 0; j < wyckoff_positions[k].size(); j++) {
        for (uint i = 0; i < wyckoff_positions[k][j].size(); i++) {
          for (uint l = 0; l < simplify(wyckoff_positions[k][j][i]).size(); l++) {
            cerr << simplify(wyckoff_positions[k][j][i]).size();
            cerr << simplify(wyckoff_positions[k][j][i])[l].dbl << simplify(wyckoff_positions[k][j][i])[l].chr << " " << endl;
          }
          cerr << wyckoff_positions[k][j][i] << " ";
        }
        cerr << endl;
      }
      cerr << endl;
    }
  }
} //namespace SYM

// ******************************************************************************
// convert_wyckoff_pos_sd
// ******************************************************************************
//Convert the wyckoff positions into sdoubles for algebraic work.
namespace SYM {
  vector<vector<vector<vector<sdouble> > > > convert_wyckoff_pos_sd(vector<vector<vector<string> > > wyckoff_positions) {
    vector<vector<vector<vector<sdouble> > > > sdw;
    //Wyckoff_positions[tuplets][singlets][singlet characters]
    for (uint k = 0; k < wyckoff_positions.size(); k++) {
      vector<vector<vector<sdouble> > > temp2;
      //cerr <<"site sizes: " <<  wyckoff_positions[k].size() << endl;
      for (uint j = 0; j < wyckoff_positions[k].size(); j++) {
        vector<vector<sdouble> > temp1;
        for (uint i = 0; i < wyckoff_positions[k][j].size(); i++) {
          //cerr << wyckoff_positions[k][j][i] << endl;
          vector<sdouble> tmpsd = simplify(wyckoff_positions[k][j][i]);
          //xb();
          //for(uint ix=0;ix<tmpsd.size();ix++){
          //  cerr <<"?: " <<  tmpsd[ix].dbl << " " << tmpsd[ix].chr << endl;
          //}
          temp1.push_back(simplify(wyckoff_positions[k][j][i]));
        }
        temp2.push_back(temp1);
      }
      //cerr << temp2.size() << endl;
      sdw.push_back(temp2);
    }
    return sdw;
  }
} //namespace SYM

// ******************************************************************************
// PointGroup_SpaceGroup (Overloaded)
// ******************************************************************************
namespace SYM {
  int* PointGroup_SpaceGroup(string pgroup) {
    int* sgrange = new int[2];
    sgrange[0] = 0;
    sgrange[1] = 0;

    if(pgroup == "1") {
      sgrange[0] = 1;
      sgrange[1] = 1;
    }
    if(pgroup == "-1") {
      sgrange[0] = 2;
      sgrange[1] = 2;
    }
    if(pgroup == "2") {
      sgrange[0] = 3;
      sgrange[1] = 5;
    }
    if(pgroup == "m") {
      sgrange[0] = 6;
      sgrange[1] = 9;
    }
    if(pgroup == "2/m") {
      sgrange[0] = 10;
      sgrange[1] = 15;
    }
    if(pgroup == "222") {
      sgrange[0] = 16;
      sgrange[1] = 24;
    }
    if(pgroup == "mm2") {
      sgrange[0] = 25;
      sgrange[1] = 46;
    }
    if(pgroup == "mmm") {
      sgrange[0] = 47;
      sgrange[1] = 74;
    }
    if(pgroup == "4") {
      sgrange[0] = 75;
      sgrange[1] = 80;
    }
    if(pgroup == "-4") {
      sgrange[0] = 81;
      sgrange[1] = 82;
    }
    if(pgroup == "4/m") {
      sgrange[0] = 83;
      sgrange[1] = 88;
    }
    if(pgroup == "422") {
      sgrange[0] = 89;
      sgrange[1] = 98;
    }
    if(pgroup == "4mm") {
      sgrange[0] = 99;
      sgrange[1] = 110;
    }
    if(pgroup == "-42m") {
      sgrange[0] = 111;
      sgrange[1] = 122;
    }
    if(pgroup == "4/mmm") {
      sgrange[0] = 123;
      sgrange[1] = 142;
    }
    if(pgroup == "3") {
      sgrange[0] = 143;
      sgrange[1] = 146;
    }
    if(pgroup == "-3") {
      sgrange[0] = 147;
      sgrange[1] = 148;
    }
    if(pgroup == "32") {
      sgrange[0] = 149;
      sgrange[1] = 155;
    }
    if(pgroup == "3m") {
      sgrange[0] = 156;
      sgrange[1] = 161;
    }
    if(pgroup == "-3m") {
      sgrange[0] = 162;
      sgrange[1] = 167;
    }
    if(pgroup == "6") {
      sgrange[0] = 168;
      sgrange[1] = 173;
    }
    if(pgroup == "-6") {
      sgrange[0] = 174;
      sgrange[1] = 174;
    }
    if(pgroup == "6/m") {
      sgrange[0] = 175;
      sgrange[1] = 176;
    }
    if(pgroup == "622") {
      sgrange[0] = 177;
      sgrange[1] = 182;
    }
    if(pgroup == "6mm") {
      sgrange[0] = 183;
      sgrange[1] = 186;
    }
    if(pgroup == "-62m" || pgroup == "-6m2") {
      sgrange[0] = 187;
      sgrange[1] = 190;
    }
    if(pgroup == "6/mmm") {
      sgrange[0] = 191;
      sgrange[1] = 194;
    }
    if(pgroup == "23") {
      sgrange[0] = 195;
      sgrange[1] = 199;
    }
    if(pgroup == "m-3") {
      sgrange[0] = 200;
      sgrange[1] = 206;
    }
    if(pgroup == "432") {
      sgrange[0] = 207;
      sgrange[1] = 214;
    }
    if(pgroup == "-43m") {
      sgrange[0] = 215;
      sgrange[1] = 220;
    }
    if(pgroup == "m-3m") {
      sgrange[0] = 221;
      sgrange[1] = 230;
    }

    return sgrange;
  }
} //namespace SYM

// ******************************************************************************
// PointGroup_SpaceGroup (Overloaded)
// ******************************************************************************
namespace SYM {
  vector<int> PointGroup_SpaceGroup(string pgroup, char c) {  //cl = centering label
    vector<int> sgs;
    int* sgrange = new int[2];
    sgrange[0] = 0;
    sgrange[1] = 0;
    //TRICLINIC
    if(pgroup == "1") {
      sgs.push_back(1);
    }
    if(pgroup == "-1") {
      sgs.push_back(2);
    }
    //MONOCLINIC
    if(pgroup == "2") {
      if(c == 'P') {
        sgs.push_back(3);
        sgs.push_back(4);
      }
      if(c == 'C') {
        sgs.push_back(5);
      }
    }
    if(pgroup == "m") {
      if(c == 'P') {
        sgs.push_back(6);
        sgs.push_back(7);
      }
      if(c == 'C') {
        sgs.push_back(8);
        sgs.push_back(9);
      }
    }
    if(pgroup == "2/m") {
      if(c == 'P') {
        sgs.push_back(10);
        sgs.push_back(11);
        sgs.push_back(13);
        sgs.push_back(14);
      }
      if(c == 'C') {
        sgs.push_back(12);
        sgs.push_back(15);
      }
    }
    //ORTHORHOMBIC
    if(pgroup == "222") {
      if(c == 'P') {
        sgs.push_back(16);
        sgs.push_back(17);
        sgs.push_back(18);
        sgs.push_back(19);
      }
      if(c == 'C') {
        sgs.push_back(20);
        sgs.push_back(21);
        sgs.push_back(21);
      }
      if(c == 'F') {
        sgs.push_back(22);
      }
      if(c == 'I') {
        sgs.push_back(23);
        sgs.push_back(24);
      }
    }
    if(pgroup == "mm2") {  //sgrange[0]=25;sgrange[1]=46;
      if(c == 'P') {
        sgs.push_back(25);
        sgs.push_back(26);
        sgs.push_back(27);
        sgs.push_back(28);
        sgs.push_back(29);
        sgs.push_back(30);
        sgs.push_back(31);
        sgs.push_back(32);
        sgs.push_back(33);
        sgs.push_back(34);
      }
      if(c == 'C') {
        sgs.push_back(35);
        sgs.push_back(36);
        sgs.push_back(37);
        sgs.push_back(38);
        sgs.push_back(39);
        sgs.push_back(40);
        sgs.push_back(41);
      }
      if(c == 'F') {
        sgs.push_back(42);
        sgs.push_back(43);
      }
      if(c == 'I') {
        sgs.push_back(44);
        sgs.push_back(45);
        sgs.push_back(46);
      }
    }

    if(pgroup == "mmm") {  //sgrange[0]=47;sgrange[1]=74;
      if(c == 'P') {
        sgs.push_back(47);
        sgs.push_back(48);
        sgs.push_back(49);
        sgs.push_back(50);
        sgs.push_back(51);
        sgs.push_back(52);
        sgs.push_back(53);
        sgs.push_back(54);
        sgs.push_back(55);
        sgs.push_back(56);
        sgs.push_back(57);
        sgs.push_back(58);
        sgs.push_back(59);
        sgs.push_back(60);
        sgs.push_back(61);
        sgs.push_back(62);
      }
      if(c == 'C') {
        sgs.push_back(63);
        sgs.push_back(64);
        sgs.push_back(65);
        sgs.push_back(66);
        sgs.push_back(67);
        sgs.push_back(68);
      }
      if(c == 'F') {
        sgs.push_back(69);
        sgs.push_back(70);
      }
      if(c == 'I') {
        sgs.push_back(71);
        sgs.push_back(72);
        sgs.push_back(73);
        sgs.push_back(74);
      }
    }
    //TETRAGONAL
    if(pgroup == "4") {  //sgrange[0]=75;sgrange[1]=80;
      if(c == 'P') {
        sgs.push_back(75);
        sgs.push_back(76);
        sgs.push_back(77);
        sgs.push_back(78);
      }
      if(c == 'I') {
        sgs.push_back(79);
        sgs.push_back(80);
      }
    }
    if(pgroup == "-4") {  //sgrange[0]=81;sgrange[1]=82;
      if(c == 'P') {
        sgs.push_back(81);
      }
      if(c == 'I') {
        sgs.push_back(82);
      }
    }
    if(pgroup == "4/m") {  //sgrange[0]=83;sgrange[1]=88;
      if(c == 'P') {
        sgs.push_back(83);
        sgs.push_back(84);
        sgs.push_back(85);
        sgs.push_back(86);
      }
      if(c == 'I') {
        sgs.push_back(87);
        sgs.push_back(88);
      }
    }
    if(pgroup == "422") {  //sgrange[0]=89;sgrange[1]=98;
      if(c == 'P') {
        sgs.push_back(89);
        sgs.push_back(90);
        sgs.push_back(91);
        sgs.push_back(92);
        sgs.push_back(93);
        sgs.push_back(94);
        sgs.push_back(95);
        sgs.push_back(96);
      }
      if(c == 'I') {
        sgs.push_back(97);
        sgs.push_back(98);
      }
    }
    if(pgroup == "4mm") {  //sgrange[0]=99;sgrange[1]=110;
      if(c == 'P') {
        sgs.push_back(99);
        sgs.push_back(100);
        sgs.push_back(101);
        sgs.push_back(102);
        sgs.push_back(103);
        sgs.push_back(104);
        sgs.push_back(105);
        sgs.push_back(106);
      }
      if(c == 'I') {
        sgs.push_back(107);
        sgs.push_back(108);
        sgs.push_back(109);
        sgs.push_back(110);
      }
    }
    if(pgroup == "-42m") {  //sgrange[0]=111;sgrange[1]=122;
      if(c == 'P') {
        sgs.push_back(111);
        sgs.push_back(112);
        sgs.push_back(113);
        sgs.push_back(114);
        sgs.push_back(115);
        sgs.push_back(116);
        sgs.push_back(117);
        sgs.push_back(118);
      }
      if(c == 'I') {
        sgs.push_back(119);
        sgs.push_back(120);
        sgs.push_back(121);
        sgs.push_back(122);
      }
    }
    if(pgroup == "4/mmm") {  //sgrange[0]=123;sgrange[1]=142;
      if(c == 'P') {
        sgs.push_back(123);
        sgs.push_back(124);
        sgs.push_back(125);
        sgs.push_back(126);
        sgs.push_back(127);
        sgs.push_back(128);
        sgs.push_back(129);
        sgs.push_back(130);
        sgs.push_back(131);
        sgs.push_back(132);
        sgs.push_back(133);
        sgs.push_back(134);
        sgs.push_back(135);
        sgs.push_back(136);
        sgs.push_back(137);
        sgs.push_back(138);
      }
      if(c == 'I') {
        sgs.push_back(139);
        sgs.push_back(140);
        sgs.push_back(141);
        sgs.push_back(142);
      }
    }
    //TRIGONAL
    if(pgroup == "3") {  //sgrange[0]=143;sgrange[1]=146;
      if(c == 'P') {
        sgs.push_back(143);
        sgs.push_back(144);
        sgs.push_back(145);
      }
      if(c == 'R') {
        sgs.push_back(146);
      }
    }
    if(pgroup == "-3") {  //sgrange[0]=147;sgrange[1]=148;
      if(c == 'P') {
        sgs.push_back(147);
      }
      if(c == 'R') {
        sgs.push_back(148);
      }
    }
    if(pgroup == "32") {  //sgrange[0]=149;sgrange[1]=155;
      if(c == 'P') {
        sgs.push_back(149);
        sgs.push_back(150);
        sgs.push_back(151);
        sgs.push_back(152);
        sgs.push_back(153);
        sgs.push_back(154);
      }
      if(c == 'R') {
        sgs.push_back(155);
      }
    }
    if(pgroup == "3m") {  //sgrange[0]=156;sgrange[1]=161;
      if(c == 'P') {
        sgs.push_back(156);
        sgs.push_back(157);
        sgs.push_back(158);
        sgs.push_back(159);
      }
      if(c == 'R') {
        sgs.push_back(160);
        sgs.push_back(161);
      }
    }
    if(pgroup == "-3m") {  //sgrange[0]=162;sgrange[1]=167;
      if(c == 'P') {
        sgs.push_back(162);
        sgs.push_back(163);
        sgs.push_back(164);
        sgs.push_back(165);
      }
      if(c == 'R') {
        sgs.push_back(166);
        sgs.push_back(167);
      }
    }
    //HEXAGONAL
    if(pgroup == "6") {  //sgrange[0]=168;sgrange[1]=173;
      sgs.push_back(168);
      sgs.push_back(169);
      sgs.push_back(170);
      sgs.push_back(171);
      sgs.push_back(172);
      sgs.push_back(173);
    }
    if(pgroup == "-6") {  //sgrange[0]=174;sgrange[1]=174;
      sgs.push_back(174);
    }
    if(pgroup == "6/m") {  //sgrange[0]=175;sgrange[1]=176;
      sgs.push_back(175);
      sgs.push_back(176);
    }
    if(pgroup == "622") {  //sgrange[0]=177;sgrange[1]=182;
      sgs.push_back(177);
      sgs.push_back(178);
      sgs.push_back(179);
      sgs.push_back(180);
      sgs.push_back(181);
      sgs.push_back(182);
    }
    if(pgroup == "6mm") {  //sgrange[0]=183;sgrange[1]=186;
      sgs.push_back(183);
      sgs.push_back(184);
      sgs.push_back(185);
      sgs.push_back(186);
    }
    if(pgroup == "-62m" || pgroup == "-6m2") {  //sgrange[0]=187;sgrange[1]=190;
      sgs.push_back(187);
      sgs.push_back(188);
      sgs.push_back(189);
      sgs.push_back(190);
    }
    if(pgroup == "6/mmm") {  //sgrange[0]=191;sgrange[1]=194;
      sgs.push_back(191);
      sgs.push_back(192);
      sgs.push_back(193);
      sgs.push_back(194);
    }
    //CUBIC
    if(pgroup == "23") {  //sgrange[0]=195;sgrange[1]=199;
      if(c == 'P') {
        sgs.push_back(195);
        sgs.push_back(198);
      }
      if(c == 'F') {
        sgs.push_back(196);
      }
      if(c == 'I') {
        sgs.push_back(197);
        sgs.push_back(199);
      }
    }
    if(pgroup == "m-3") {  //sgrange[0]=200;sgrange[1]=206;
      if(c == 'P') {
        sgs.push_back(200);
        sgs.push_back(201);
        sgs.push_back(205);
      }
      if(c == 'F') {
        sgs.push_back(202);
        sgs.push_back(203);
      }
      if(c == 'I') {
        sgs.push_back(204);
        sgs.push_back(206);
      }
    }
    if(pgroup == "432") {  //sgrange[0]=207;sgrange[1]=214;
      if(c == 'P') {
        sgs.push_back(207);
        sgs.push_back(208);
        sgs.push_back(212);
        sgs.push_back(213);
      }
      if(c == 'F') {
        sgs.push_back(209);
        sgs.push_back(210);
      }
      if(c == 'I') {
        sgs.push_back(211);
        sgs.push_back(214);
      }
    }
    if(pgroup == "-43m") {  //sgrange[0]=215;sgrange[1]=220;
      if(c == 'P') {
        sgs.push_back(215);
        sgs.push_back(218);
      }
      if(c == 'F') {
        sgs.push_back(216);
        sgs.push_back(219);
      }
      if(c == 'I') {
        sgs.push_back(217);
        sgs.push_back(220);
      }
    }
    if(pgroup == "m-3m") {  //sgrange[0]=221;sgrange[1]=230;
      if(c == 'P') {
        sgs.push_back(221);
        sgs.push_back(222);
        sgs.push_back(223);
        sgs.push_back(224);
      }
      if(c == 'F') {
        sgs.push_back(225);
        sgs.push_back(226);
        sgs.push_back(227);
        sgs.push_back(228);
      }
      if(c == 'I') {
        sgs.push_back(229);
        sgs.push_back(230);
      }
    }

    return sgs;
  }
} //namespace SYM

// **********************************************************************************************************************
// getLatticeVectorsFromOriginalMirrorOperations
// **********************************************************************************************************************
// Obtain new lattice vectors corresponding to mirrors
namespace SYM {
  vector<xvector<double> > getLatticeVectorsFromOriginalMirrorOperations(vector<Glide>& old_mirrors, vector<Glide>& new_mirrors,
      vector<xvector<double> >& lattice_vectors,
      bool& all_matched) {
    xvector<double> null_tol;
    null_tol(1) = _ZERO_TOL_;
    null_tol(2) = _ZERO_TOL_;
    null_tol(3) = _ZERO_TOL_;
    vector<xvector<double> > new_lattice_vectors;
    for (uint i = 0; i < new_mirrors.size(); i++) {
      bool found = false;
      for (uint j = 0; j < old_mirrors.size(); j++) {
        xvector<double> dir_diff = new_mirrors[i].return_direction() - old_mirrors[j].return_direction();
        xvector<double> dir_add = new_mirrors[i].return_direction() + old_mirrors[j].return_direction();
        xvector<double> point_diff = new_mirrors[i].return_point() - old_mirrors[j].return_point();
        bool diff = (aurostd::abs(dir_diff(1)) < _ZERO_TOL_ && aurostd::abs(dir_diff(2)) < _ZERO_TOL_ && aurostd::abs(dir_diff(3)) < _ZERO_TOL_);
        bool add = (aurostd::abs(dir_add(1)) < _ZERO_TOL_ && aurostd::abs(dir_add(2)) < _ZERO_TOL_ && aurostd::abs(dir_add(3)) < _ZERO_TOL_);
        bool point = (aurostd::abs(point_diff(1)) < _ZERO_TOL_ && aurostd::abs(point_diff(2)) < _ZERO_TOL_ && aurostd::abs(point_diff(3)) < _ZERO_TOL_);
        if((diff || add) && point) {
          if(add) {
            new_lattice_vectors.push_back(-lattice_vectors[j]);
          } else {
            new_lattice_vectors.push_back(lattice_vectors[j]);
          }
          found = true;
          break;
        }
      }
      if(found == false) {
        all_matched = false;
      }
    }
    return new_lattice_vectors;
  }
} //namespace SYM 

// **********************************************************************************************************************
// getLatticeVectorsFromOriginalRotationOperations
// **********************************************************************************************************************
// Obtain new lattice vectors corresponding to rotations
namespace SYM {
  vector<xvector<double> > getLatticeVectorsFromOriginalRotationOperations(vector<Screw>& old_rotations_twofold,
      vector<Screw>& old_rotations_rot, vector<Screw>& new_rotations,
      vector<xvector<double> >& twofold_lattice_vectors,
      vector<xvector<double> >& rot_lattice_vectors,
      bool& all_matched) {
    vector<xvector<double> > new_lattice_vectors;
    for (uint i = 0; i < new_rotations.size(); i++) {
      bool found = false;
      for (uint j = 0; j < old_rotations_twofold.size(); j++) {
        xvector<double> dir_diff = new_rotations[i].return_direction() - old_rotations_twofold[j].return_direction();
        xvector<double> dir_add = new_rotations[i].return_direction() + old_rotations_twofold[j].return_direction();
        xvector<double> point_diff = new_rotations[i].return_point() - old_rotations_twofold[j].return_point();
        bool diff = (aurostd::abs(dir_diff(1)) < _ZERO_TOL_ && aurostd::abs(dir_diff(2)) < _ZERO_TOL_ && aurostd::abs(dir_diff(3)) < _ZERO_TOL_);
        bool add = (aurostd::abs(dir_add(1)) < _ZERO_TOL_ && aurostd::abs(dir_add(2)) < _ZERO_TOL_ && aurostd::abs(dir_add(3)) < _ZERO_TOL_);
        bool point = (aurostd::abs(point_diff(1)) < _ZERO_TOL_ && aurostd::abs(point_diff(2)) < _ZERO_TOL_ && aurostd::abs(point_diff(3)) < _ZERO_TOL_);
        if((diff || add) && point) {
          if(add) {
            new_lattice_vectors.push_back(-twofold_lattice_vectors[j]);
          } else {
            new_lattice_vectors.push_back(twofold_lattice_vectors[j]);
          }
          found = true;
          break;
        }
      }
      if(found == false) {
        for (uint j = 0; j < old_rotations_rot.size(); j++) {
          xvector<double> dir_diff = new_rotations[i].return_direction() - old_rotations_rot[j].return_direction();
          xvector<double> dir_add = new_rotations[i].return_direction() + old_rotations_rot[j].return_direction();
          xvector<double> point_diff = new_rotations[i].return_point() - old_rotations_rot[j].return_point();
          bool diff = (aurostd::abs(dir_diff(1)) < _ZERO_TOL_ && aurostd::abs(dir_diff(2)) < _ZERO_TOL_ && aurostd::abs(dir_diff(3)) < _ZERO_TOL_);
          bool add = (aurostd::abs(dir_add(1)) < _ZERO_TOL_ && aurostd::abs(dir_add(2)) < _ZERO_TOL_ && aurostd::abs(dir_add(3)) < _ZERO_TOL_);
          bool point = (aurostd::abs(point_diff(1)) < _ZERO_TOL_ && aurostd::abs(point_diff(2)) < _ZERO_TOL_ && aurostd::abs(point_diff(3)) < _ZERO_TOL_);
          if((diff || add) && point) {
            if(add) {
              new_lattice_vectors.push_back(-rot_lattice_vectors[j]);
            } else {
              new_lattice_vectors.push_back(rot_lattice_vectors[j]);
            }
            found = true;
            break;
          }
        }
      }
      if(found == false) {
        all_matched = false;
      }
    }
    return new_lattice_vectors;
  }
} //namespace SYM

// ******************************************************************************
// Fold atoms into cell
// ******************************************************************************
// Folds atoms into the cell
//CO20180409 - moved to xatom
//[CO20180409 - moved to xatom]namespace SYM {
//[CO20180409 - moved to xatom]  deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool& skew) {
//[CO20180409 - moved to xatom]    double tol = _SYM_TOL_;
//[CO20180409 - moved to xatom]    return foldAtomsInCell(atoms, c2f_new, f2c_new, skew, tol);
//[CO20180409 - moved to xatom]  }
//[CO20180409 - moved to xatom]} //namespace SYM

//[CO20180409 - moved to xatom]namespace SYM {
//[CO20180409 - moved to xatom]  deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool& skew, double& tol) {
//[CO20180409 - moved to xatom]    deque<_atom> atoms_in_cell;
//[CO20180409 - moved to xatom]    _atom tmp;
//[CO20180409 - moved to xatom]    for (uint j = 0; j < atoms.size(); j++) {
//[CO20180409 - moved to xatom]      if(atoms_in_cell.size() == 0) {
//[CO20180409 - moved to xatom]	//[OBSOLETE]atoms[j].fpos = c2f_new * atoms[j].cpos;
//[CO20180409 - moved to xatom]	//[OBSOLETE]atoms[j].fpos = BringInCell(atoms[j].fpos);
//[CO20180409 - moved to xatom]	//[OBSOLETE]atoms[j].cpos = f2c_new * atoms[j].fpos;
//[CO20180409 - moved to xatom]	atoms_in_cell.push_back(atoms[j]);
//[CO20180409 - moved to xatom]	atoms_in_cell.back().fpos = BringInCell(c2f_new * atoms[j].cpos);
//[CO20180409 - moved to xatom]  atoms_in_cell.back().cpos = f2c_new * atoms_in_cell.back().fpos;
//[CO20180409 - moved to xatom]  atoms_in_cell.back().ijk(1)=0; atoms_in_cell.back().ijk(2)=0; atoms_in_cell.back().ijk(3)=0;
//[CO20180409 - moved to xatom]      } else {
//[CO20180409 - moved to xatom]  //bool duplicate_atom = false;
//[CO20180409 - moved to xatom]  tmp.fpos = BringInCell(c2f_new * atoms[j].cpos);
//[CO20180409 - moved to xatom]  tmp.cpos = f2c_new * tmp.fpos;
//[CO20180409 - moved to xatom]  //[OBSOLETE]for (uint a = 0; a < atoms_in_cell.size(); a++) {
//[CO20180409 - moved to xatom]  //[OBSOLETE]  if(MapAtomsInNewCell(atoms_in_cell[a], tmp, c2f_orig, f2c_new, skew, tol)) {
//[CO20180409 - moved to xatom]  //[OBSOLETE]  if(MapAtoms(atoms_in_cell[a], tmp, c2f_orig, f2c_new, skew, tol)) {
//[CO20180409 - moved to xatom]  //[OBSOLETE]    duplicate_atom = true;
//[CO20180409 - moved to xatom]  //[OBSOLETE]    break;
//[CO20180409 - moved to xatom]  //[OBSOLETE]  }
//[CO20180409 - moved to xatom]  //[OBSOLETE]}
//[CO20180409 - moved to xatom]  //[OBSOLETE]if(duplicate_atom == false) {
//[CO20180409 - moved to xatom]  if(!SYM::MapAtom(atoms_in_cell,tmp,false,c2f_new,f2c_new,skew,tol)){
//[CO20180409 - moved to xatom]    //[OBSOLETE]atoms[j].fpos = tmp.fpos; //BringInCell(tmp.fpos);
//[CO20180409 - moved to xatom]    //[OBSOLETE]atoms[j].cpos = tmp.cpos; //f2c_new * atoms[j].fpos;
//[CO20180409 - moved to xatom]    atoms_in_cell.push_back(atoms[j]);
//[CO20180409 - moved to xatom]    atoms_in_cell.back().fpos = tmp.fpos;
//[CO20180409 - moved to xatom]    atoms_in_cell.back().cpos = tmp.cpos;
//[CO20180409 - moved to xatom]    atoms_in_cell.back().ijk(1)=0; atoms_in_cell.back().ijk(2)=0; atoms_in_cell.back().ijk(3)=0;
//[CO20180409 - moved to xatom]  }
//[CO20180409 - moved to xatom]      }
//[CO20180409 - moved to xatom]    }
//[CO20180409 - moved to xatom]    return atoms_in_cell;
//[CO20180409 - moved to xatom]  }
//[CO20180409 - moved to xatom]} //namespace SYM

// ******************************************************************************
// Map Atoms into a New Cell
// ******************************************************************************
// This differs from MapAtom since we need information about the original and new lattice
namespace SYM {
  bool MapAtomsInNewCell(_atom& a, _atom& b, xmatrix<double>& lattice_new, bool& skew, double& tol) { //DX20190619 - changed c2f_orig and f2c_new to lattice_new
    //DX20190613 [OBSOLETE-moved down] xvector<double> fdiff = a.fpos - b.fpos;
    double min_dist = 1e9;
    if(skew) {
      xvector<double> min_cdiff = SYM::minimizeDistanceCartesianMethod(a.cpos, b.cpos, lattice_new); //DX20190613
      min_dist = aurostd::modulus(min_cdiff); //DX20190613
      //DX20190613 [OBSOLETE] return SYM::minimizeCartesianDistance(a.cpos, b.cpos, fdiff, c2f_orig, f2c_new, tol);
    } else {
      xmatrix<double> f2c_new = trasp(lattice_new); //DX20190613
      xvector<double> fdiff = a.fpos - b.fpos; //DX20190613
      xvector<double> min_fdiff = SYM::minimizeDistanceFractionalMethod(fdiff); //DX20190613 
      min_dist = aurostd::modulus(f2c_new * min_fdiff); //DX20190613
      //DX20190613 [OBSOLETE] SYM::PBC(fdiff);
    }
    return (min_dist < tol); //DX20190613
    //DX20190613 [OBSOLETE] return (aurostd::modulus(f2c_new * fdiff) < tol);
  }
} //namespace SYM

// ******************************************************************************
// Map Atoms into a New Cell
// ******************************************************************************
// This differs from MapAtom since we need information about the original and new lattice
namespace SYM {
  bool MapAtomsInNewCell(xvector<double>& a, xvector<double>& b, xmatrix<double>& lattice_new, bool& skew, double& tol) { //DX20190619 - changed c2f_orig and f2c_new to lattice_new
    //DX20190613 [OBSOLETE-moved down] xvector<double> fdiff = a - b;
    double min_dist = 1e9; //DX20190613
    xmatrix<double> f2c_new = trasp(lattice_new); //DX20190613
    if(skew) {
      xvector<double> cpos1 = f2c_new * a; //DX20190613 - coord1 to cpos1
      xvector<double> cpos2 = f2c_new * b; //DX20190613 - coord2 to cpos2
      xvector<double> min_cdiff = SYM::minimizeDistanceCartesianMethod(cpos1, cpos2, lattice_new); //DX20190613
      min_dist = aurostd::modulus(min_cdiff);
      //DX20190613 [OBSOLETE] return SYM::minimizeCartesianDistance(coord1, coord2, fdiff, c2f_orig, f2c_new, tol);
    } else {
      xvector<double> fdiff = a - b; //DX20190613
      xvector<double> min_fdiff = SYM::minimizeDistanceFractionalMethod(fdiff); //DX20190613
      min_dist = aurostd::modulus(f2c_new * min_fdiff); //DX20190613
      //DX20190613 [OBSOLETE] SYM::PBC(fdiff);
    }
    return (min_dist < tol); //DX20190613
    //DX20190613 [OBSOLETE] return (aurostd::modulus(f2c_new * fdiff) < tol);
  }
} //namespace SYM

// ******************************************************************************
// Group symmetry equivalent atoms
// ******************************************************************************
// Group symmetrically equivalent atoms
namespace SYM {
  deque<deque<_atom> > groupSymmetryEquivalentAtoms(deque<_atom>& atoms, xmatrix<double>& lattice, vector<xmatrix<double> >& sym_ops,
      vector<xvector<double> >& translations, double& min_dist, double& tol) { //DX20190215 - added tol
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    deque<deque<_atom> > equivalent_atoms;
    xmatrix<double> f2c = trasp(lattice);
    //DX20190619 [OBSOLETE] xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = SYM::isLatticeSkewed(lattice, min_dist, tol); //DX20190215
    for (uint k = 0; k < atoms.size(); k++) {
      deque<_atom> tmpv;
      for (uint i = 0; i < sym_ops.size(); i++) {
        _atom tmp = atoms[k]; //DX20191011 - initialized instead of setting individually below
        //DX20190905 [OBSOLETE-no more mod_one_xvec] tmp.fpos = (SYM::mod_one_xvec((sym_ops[i] * atoms[k].fpos + translations[i])));
        tmp.fpos = sym_ops[i] * atoms[k].fpos + translations[i]; //DX20190905
        BringInCellInPlace(tmp.fpos); //DX20190905
        tmp.cpos = f2c * tmp.fpos;
        //DX20191011 [OBSOLETE] tmp.name = atoms[k].name;
        //DX20191011 [OBSOLETE] tmp.type = atoms[k].type;
        //DX20191011 [OBSOLETE] tmp.spin = atoms[k].spin; //DX20170921 - magnetic sym
        //DX20191011 [OBSOLETE] tmp.spin_is_given = atoms[k].spin_is_given; //DX20170921 - magnetic sym
        //DX20191011 [OBSOLETE] tmp.noncoll_spin = atoms[k].noncoll_spin; //DX20171205 - magnetic sym (non-collinear)
        //DX20191011 [OBSOLETE] tmp.noncoll_spin_is_given = atoms[k].noncoll_spin_is_given; //DX20171205 - magnetic sym (non-collinear)
        bool contained = false;
        for (uint j = 0; j < equivalent_atoms.size(); j++) {
          if(SYM::MapAtom(equivalent_atoms[j], tmp, TRUE, lattice, f2c, skew, tol)) { //DX20190215 //DX20190619 - lattice and f2c as input
            contained = true;
            break;
          }
        }
        if(contained == false) {
          tmpv.push_back(tmp);
        }
      }
      if(tmpv.size() > 0) {
        reduce_atom_deques(tmpv, lattice, min_dist,tol); //DX20190215
        equivalent_atoms.push_back(tmpv);
      }
    }
    // ===== DEBUG::PRINT EQUIVALENT ATOMS ===== //
    if(LDEBUG) {
      cerr << "SYM::shiftSymmetryEquivalentAtoms: Equivalent atoms: " << endl;
      for(uint i=0;i<equivalent_atoms.size();i++){
        xb();
        print(equivalent_atoms[i]);
      }
      cerr << endl;
    }
    return equivalent_atoms;
  }
} //namespace SYM

// ******************************************************************************
// Shift symmetry equivalent atoms
// ******************************************************************************
// Shift symmetrically equivalent atoms by an origin shift
namespace SYM {
  deque<deque<_atom> > shiftSymmetryEquivalentAtoms(deque<deque<_atom> >& equivalent_atoms, xmatrix<double>& lattice, xvector<double>& translation, double& min_dist, double& tol) { //DX20190215
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    deque<deque<_atom> > equivalent_atoms_shifted;
    xmatrix<double> f2c = trasp(lattice);
    //DX20190619 [OBSOLETE] xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = SYM::isLatticeSkewed(lattice, min_dist, tol); //DX20190215

    deque<_atom> one_shifted_group;
    for (uint i = 0; i < equivalent_atoms.size(); i++) {
      one_shifted_group.clear();
      for (uint j = 0; j < equivalent_atoms[i].size(); j++) {
        _atom tmp = equivalent_atoms[i][j]; //DX20191011 - initialized instead of setting individually below
        //DX20190905 [OBSOLETE-no more mod_one_xvec] xvector<double> tmpxvec = SYM::mod_one_xvec(equivalent_atoms[i][j].fpos + translation);
        xvector<double> tmpxvec = equivalent_atoms[i][j].fpos + translation; //DX20190905
        BringInCellInPlace(tmpxvec); //DX20190905
        // ===== Ensure that multiple equivalent atoms are not occupying the same space ===== //
        tmp.fpos = tmpxvec;
        tmp.cpos = f2c * tmp.fpos;
        //DX20191011 [OBSOLETE] tmp.name = equivalent_atoms[i][j].name;
        //DX20191011 [OBSOLETE] tmp.type = equivalent_atoms[i][j].type;
        //DX20191011 [OBSOLETE] tmp.spin = equivalent_atoms[i][j].spin; //DX20170921 - magnetic sym
        //DX20191011 [OBSOLETE] tmp.spin_is_given = equivalent_atoms[i][j].spin_is_given; //DX20170921 - magnetic sym
        //DX20191011 [OBSOLETE] tmp.noncoll_spin = equivalent_atoms[i][j].noncoll_spin; //DX20171205 - magnetic sym (non-collinear)
        //DX20191011 [OBSOLETE] tmp.noncoll_spin_is_given = equivalent_atoms[i][j].noncoll_spin_is_given; //DX20171205 - magnetic sym (non-collinear)
        if(one_shifted_group.size() == 0) {
          one_shifted_group.push_back(tmp);
        } else {
          bool duplicate_atom = false;
          for (uint b = 0; b < one_shifted_group.size(); b++) {
            if(SYM::MapAtom(one_shifted_group[b].fpos, tmp.fpos, lattice, f2c, skew, tol)) { //DX20190215 //DX20190619 - lattice and f2c as input
              duplicate_atom = true;
              break;
            }
          }
          if(duplicate_atom == false) {
            one_shifted_group.push_back(tmp);
          }
        }
      }
      equivalent_atoms_shifted.push_back(one_shifted_group);
    }
    // ===== DEBUG::PRINT EQUIVALENT ATOMS ===== //
    if(LDEBUG) {
      cerr << "SYM::shiftSymmetryEquivalentAtoms: Equivalent atoms after origin shift (" << translation << "): " << endl;
      for(uint i=0;i<equivalent_atoms_shifted.size();i++){
        xb();
        print(equivalent_atoms_shifted[i]);
      }
      cerr << endl;
    }
    return equivalent_atoms_shifted;
  }
} //namespace SYM

// ******************************************************************************
// GCD_conventional_atomic_basis
// ******************************************************************************
//Determine if the conventional atomic basis ratio is consistent with the primitive atomic basis
namespace SYM {
  bool GCD_conventional_atomic_basis(deque<_atom>& conventional_basis_atoms, deque<deque<_atom> >& prim_split_atom_types, int& prim_GCD) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    deque<deque<_atom> > split_atom_types = SYM::break_up_by_type(conventional_basis_atoms);
    if(split_atom_types.size() != prim_split_atom_types.size()) {
      if(LDEBUG) {
        cerr << "SYM::GCD_conventional_atomic_basis::WARNING: Number of atom types is not the same: " << split_atom_types.size() << " != " << prim_split_atom_types.size() << endl;
      }
      return false;
    }
    bool consistent_ratio = true;
    if(split_atom_types.size()==0){throw aurostd::xerror(__AFLOW_FILE__,"SYM::GCD_conventional_atomic_basis():","split_atom_types.size()==0",_INPUT_ERROR_);}
    int GCD_num = 1;  //CO+DX20191201
    if(split_atom_types.size() > 1) { // If only one type of atom, set to 1
      GCD_num=(int)split_atom_types[0].size();  //CO+DX20191201
      for (uint p = 0; p < split_atom_types.size() - 1; p++) {
        aurostd::GCD(GCD_num, (int)split_atom_types[p + 1].size(),GCD_num); //CO20191201
      }
      for (uint p = 0; p < prim_split_atom_types.size(); p++) {
        for (uint s = 0; s < split_atom_types.size(); s++) {
          if(prim_split_atom_types[p][0].name == split_atom_types[s][0].name) {
            if((prim_split_atom_types[p].size() * GCD_num) / prim_GCD != split_atom_types[s].size()) {
              consistent_ratio = false;
            }
          }
        }
      }
    } else if(split_atom_types.size() == 0) {  //Conventional basis atoms == 0
      if(LDEBUG) { cerr << "SYM::GCD_conventional_atomic_basis::WARNING: Conventional Basis atoms = 0" << endl; }
      consistent_ratio = false;
    }
    //cerr << "consistent_ratio: " << consistent_ratio << endl;
    return consistent_ratio;
  }
} //namespace SYM

// ******************************************************************************
// invecindex
// ******************************************************************************
//Vector search functions
namespace SYM {
  int invecindex(vector<_atom> vec, _atom n) {
    int out = 0;
    for (uint i = 0; i < vec.size(); i++) {
      if(vec[i].fpos == n.fpos && vec[i].name == n.name) {
        out = i;
      }
    }
    return out;
  }
} //namespace SYM

//TOPOLOGY FUNCTIONS
// **********************************************************************************************************************
// find_vectors_inplane
// **********************************************************************************************************************
// Finds the vectors perpendicular to a given vector (i.e. in a given plane)
namespace SYM {
  vector<xvector<double> > find_vectors_inplane(const vector<xvector<double> >& big_expanded, const xvector<double>& perp_to_vec, double& tol) { //DX20190215
    //DX20190215 [OBSOLETE] double tol = _SYM_TOL_;
    vector<xvector<double> > vectors_in_plane;
    vector<double> pnorms;
    vector<int> indicies_vecs_in_plane;
    for (uint i = 0; i < big_expanded.size(); i++) {
      if(aurostd::modulus(big_expanded[i]) > _ZERO_TOL_ && aurostd::modulus(perp_to_vec) > _ZERO_TOL_) {
        if(aurostd::abs(((aurostd::modulus(big_expanded[i]) + aurostd::modulus(perp_to_vec)) / (2.0)) * sin(aurostd::angle(big_expanded[i], perp_to_vec) - (Pi_r / 2.0))) < tol) {
          pnorms.push_back(aurostd::modulus(big_expanded[i]));
          indicies_vecs_in_plane.push_back(i);
        }
      }
    }
    vector<int> min_vecs;
    for (uint i = 0; i < indicies_vecs_in_plane.size(); i++) {
      min_vecs.push_back(i);
    }
    for (uint i = 0; i < min_vecs.size(); i++) {
      vectors_in_plane.push_back(big_expanded[indicies_vecs_in_plane[min_vecs[i]]]);
    }
    return vectors_in_plane;
  }
} //namespace SYM

// **********************************************************************************************************************
// count_types
// **********************************************************************************************************************
//counts the types of atoms and outputs vector<int>
namespace SYM {
  vector<int> count_types(deque<_atom>& vatom) {
    vector<int> types_count;
    deque<deque<_atom> > vaa = SYM::break_up_by_type(vatom);
    for (uint i = 0; i < vaa.size(); i++) {
      types_count.push_back(vaa[i].size());
    }
    return types_count;
  }
} //namespace SYM

// **********************************************************************************************************************
// SYM::countWyckoffTypes()
// **********************************************************************************************************************
//counts the types of Wyckoff positions per species
namespace SYM {
  vector<int> countWyckoffTypes(const vector<wyckoffsite_ITC>& Wyckoff_positions){
    vector<int> types, types_count;
    size_t index = 0;
    for(size_t i=0;i<Wyckoff_positions.size();i++){
      if(!aurostd::WithinList(types,(int)Wyckoff_positions[i].index,index)){
        types.push_back(Wyckoff_positions[i].index);
        types_count.push_back(1);
      }
      else{ types_count[index] +=1; }
    }
    return types_count;
  }
} //namespace SYM

// **********************************************************************************************************************
// mirror_operations
// **********************************************************************************************************************
//GET MIRROR OPERATIONS
namespace SYM {
  vector<Glide> mirror_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew, double& tol) { //DX20190215
    //GET MIRROR PLANES
    //*********************************************************************
    // This part is a little tricky because planes that are symmetrically unique are
    // excluded due to the way symmetry_axes_equivalent works: For example, the planes going
    // through 0,0,0 and 0,0,.5 with the same normal will be excluded, even though these
    // are not lattice-translates. However, this doesn't seem to be a problem because of the
    // way the point groups are defined group-theoretically. Keep in mind, though, that when I
    // get the triplet-based operations, these pairs are not eliminated, resulting in more
    // operations than are usually listed
    //*********************************************************************

    vector<Glide> mirror_ops_vec;
    uint tsize = expanded_lattice.size();
    uint doubles_size = CombinationNr(2, tsize);
    vector<int> cmb_dbl = AllCombination41(2, tsize, 1);
    vector<xvector<double> > tmp_dbl;
    xvector<double> mpnorm;
    xvector<double> mpmid;
    xvector<double> nullvec;
    for (uint i = 0; i < doubles_size; i++) {
      tmp_dbl.clear();
      tmp_dbl.push_back(expanded_lattice[cmb_dbl[0]]);
      tmp_dbl.push_back(expanded_lattice[cmb_dbl[1]]);
      Glide candidate_mirror;
      mpnorm = tmp_dbl[1] - tmp_dbl[0];
      mpmid = tmp_dbl[1] + (1.0 / 2.0) * mpnorm;
      mpnorm = (1.0 / aurostd::modulus(mpnorm)) * mpnorm;
      candidate_mirror.get_glide_direct(mpnorm, nullvec);  //mpmid);
      bool cont_dbl = true;
      xvector<double> lattice_vector = tmp_dbl[1] - tmp_dbl[0];
      for (int i = 0; i < 8; i++) {
        //if(!is_lattice_point(Linv,candidate_mirror*expanded_cell[i]))
        if(!is_lattice_point(L, candidate_mirror * expanded_cell[i], lattice_vector, radius, skew, tol)) //DX20190215 - added tol
        { //CO20200106 - patching for auto-indenting
          cont_dbl = false;
        }
      }
      if(cont_dbl == true) {
        is_lattice_point(L, tmp_dbl[1] - tmp_dbl[0], lattice_vector, radius, skew, tol); //DX20190215 - added tol
        bool contained_in_list = false;
        for (uint i = 0; i < mirror_ops_vec.size(); i++) {
          //Use function to compare rotation axes
          if(symmetry_axis_equivalent(L, Linv, nullvec, mpnorm, mirror_ops_vec[i].return_point(), mirror_ops_vec[i].return_direction(), tol)) { //DX20190215 - added tol
            contained_in_list = true;
            if(aurostd::modulus(lattice_vector) < aurostd::modulus(lattice_vectors[i])) {
              lattice_vectors[i] = lattice_vector;
            }
            break;
          }
        }
        //Add to list if not already there
        if(contained_in_list == false) {
          mirror_ops_vec.push_back(candidate_mirror);
          lattice_vectors.push_back(lattice_vector);
        }
      }
      cmb_dbl = AllCombination42(2, tsize, cmb_dbl);
    }
    return mirror_ops_vec;
  }
} //namespace SYM

// **********************************************************************************************************************
// twofold_operations
// **********************************************************************************************************************
namespace SYM {
  vector<Screw> twofold_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew, double& tol) { //DX20190215 - added tol
    xmatrix<double> f2c = trasp(L);
    xmatrix<double> c2f = inverse(trasp(L));

    vector<Screw> twofold_ops_vec;
    uint tsize = expanded_lattice.size();
    const xvector<double> Nullvec;
    int doubles_size = CombinationNr(2, tsize);
    vector<int> cmb_dbl = AllCombination41(2, tsize, 1);
    vector<xvector<double> > tmp_dbl;
    for (int i = 0; i < doubles_size; i++) {
      tmp_dbl.clear();
      tmp_dbl.push_back(expanded_lattice[cmb_dbl[0]]);
      tmp_dbl.push_back(expanded_lattice[cmb_dbl[1]]);
      Screw candidate_rotation;
      xvector<double> axis_direction = tmp_dbl[1] - tmp_dbl[0];
      axis_direction = (1 / aurostd::modulus(axis_direction)) * axis_direction;
      xvector<double> point_on_axis = Nullvec;
      candidate_rotation.get_screw_direct(axis_direction, point_on_axis, 2);
      bool cont = true;
      xvector<double> lattice_vector;
      for (int i = 0; i < 8; i++) {
        if(!points_equivalent(c2f, f2c, Nullvec, candidate_rotation * expanded_cell[i], lattice_vector, radius, skew, tol)) { //DX20190215 - added tol
          cont = false;
        }
      }
      if(cont == true) {
        points_equivalent(c2f, f2c, Nullvec, tmp_dbl[1] - tmp_dbl[0], lattice_vector, radius, skew, tol); //DX20190215
        //MUST ELIMINATE DUPLICATES THAT ARISE BECAUSE OF EXPANDED LATTICE
        //Check if operation is included in list
        bool contained_in_list = false;
        for (uint i = 0; i < twofold_ops_vec.size(); i++) {
          //Use function to compare rotation axes
          if(symmetry_axis_equivalent(L, Linv, point_on_axis, axis_direction, twofold_ops_vec[i].return_point(), twofold_ops_vec[i].return_direction(), tol)) { //DX20190215 - added tol
            if(aurostd::modulus(lattice_vector) < aurostd::modulus(lattice_vectors[i])) {
              lattice_vectors[i] = lattice_vector;
            }
            contained_in_list = true;
            break;
          }
        }
        //Add to list if not already there
        if(contained_in_list == false) {
          twofold_ops_vec.push_back(candidate_rotation);
          lattice_vectors.push_back(lattice_vector);
        }
      }
      cmb_dbl = AllCombination42(2, tsize, cmb_dbl);
    }
    return twofold_ops_vec;
  }
} //namespace SYM

// **********************************************************************************************************************
// triplet_operations
// **********************************************************************************************************************
namespace SYM {
  vector<Screw> triplet_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew, double& tol) { //DX20190215
    int q = 3;
    vector<xvector<double> > big_expanded = expand_lattice_positive_only(q, q, q, L);

    xmatrix<double> f2c = trasp(L);
    xmatrix<double> c2f = inverse(trasp(L));

    vector<Screw> rotation_ops_vec;
    double rotation_order = 1.0; //DX20180514 - added initialization
    uint tsize = expanded_lattice.size();
    xvector<double> Nullvec;

    int triplets_size = CombinationNr(3, tsize);
    vector<int> cmb = AllCombination41(3, tsize, 1);
    vector<xvector<double> > tmp_trip;
    vector<vector<xvector<double> > > point_sets;
    char tmpchar;
    int count = 0;
    int rotation_ops = 0;
    for (int i = 0; i < triplets_size; i++) {
      tmp_trip.clear();
      tmp_trip.push_back(expanded_lattice[cmb[0]]);
      tmp_trip.push_back(expanded_lattice[cmb[1]]);
      tmp_trip.push_back(expanded_lattice[cmb[2]]);
      tmpchar = discern_rot_sym(get_mod_angle_tensor(tmp_trip), tol); //DX20190215 - added tol
      if(tmpchar != 'z') {
        count++;
        Screw candidate_rotation;
        //Get axis direction
        xvector<double> axis_direction = aurostd::vector_product(tmp_trip[1] - tmp_trip[0], tmp_trip[2] - tmp_trip[0]); //DX20190905 - SYM::CrossPro() -> aurostd::vector_product()
        xvector<double> full_direction = axis_direction;
        bool found_point = false;
        for (uint j = 0; j < big_expanded.size(); j++) {
          for (uint k = j + 1; k < big_expanded.size(); k++) {
            if((((aurostd::modulus(axis_direction) + aurostd::modulus(big_expanded[j] - big_expanded[k])) / 2.0) * sin(aurostd::angle(axis_direction, big_expanded[j] - big_expanded[k])) < tol)) { //DX20190215
              full_direction = big_expanded[j] - big_expanded[k];
              axis_direction = big_expanded[j] - big_expanded[k];
              found_point = true;
              break;
            }
          }
          if(found_point == true) {
            break;
          }
        }
        if(found_point == false) {
          cmb = AllCombination42(3, tsize, cmb);
          continue;
        }
        axis_direction = (1 / aurostd::modulus(axis_direction)) * axis_direction;
        //The point on axis can always be 0,0,0 for lattice (equivalent points and group theory)
        xvector<double> point_on_axis = Nullvec;
        //double zero_one = aurostd::modulus(tmp_trip[0]-tmp_trip[1]);
        //double zero_two = aurostd::modulus(tmp_trip[0]-tmp_trip[2]);
        if(tmpchar == 'a') {  //3-fold
          //point_on_axis = (1.0/3.0)*(tmp_trip[0]+tmp_trip[1]+tmp_trip[2]);
          candidate_rotation.get_screw_direct(axis_direction, point_on_axis, 3);
          rotation_order = 3.0;
        }
        if(tmpchar == 'b') {  //4-fold
          //if((zero_one-zero_two)>tol){
          //point_on_axis = (.5)*(tmp_trip[0]+tmp_trip[1]);
          //}
          //if((zero_two-zero_one)>tol){
          //point_on_axis = (.5)*(tmp_trip[0]+tmp_trip[2]);
          //}
          //if(aurostd::abs(zero_one-zero_two)<tol){
          //point_on_axis = (.5)*(tmp_trip[1]+tmp_trip[2]);
          //}
          candidate_rotation.get_screw_direct(axis_direction, point_on_axis, 4);
          rotation_order = 4.0;
        }
        if(tmpchar == 'c') {  //6-fold
          //if((zero_one-zero_two)>tol){
          //point_on_axis = (tmp_trip[0]+tmp_trip[1])-tmp_trip[2];
          //}
          //if((zero_two-zero_one)>tol){
          //point_on_axis = (tmp_trip[0]+tmp_trip[2])-tmp_trip[1];
          //}
          //if(aurostd::abs(zero_one-zero_two)<tol){
          //point_on_axis = (tmp_trip[1]+tmp_trip[2])-tmp_trip[0];
          //}
          candidate_rotation.get_screw_direct(axis_direction, point_on_axis, 6);
          rotation_order = 6.0;
        }
        bool cont = true;
        xvector<double> lattice_vector;
        //cerr << "candidate_rotation: " << candidate_rotation << endl;
        //cerr << "+++++++++++++++++++++++++++++++++++++++" << endl;
        for (int i = 0; i < 8; i++) {
          if(!points_equivalent(c2f, f2c, Nullvec, candidate_rotation * expanded_cell[i], lattice_vector, radius, skew, tol)) { //DX20190215 - added tol
            cont = false;
          }
        }
        if(cont == true) {
          points_equivalent(c2f, f2c, Nullvec, full_direction, lattice_vector, radius, skew, tol); //DX20190215 - added tol
          //MUST ELIMINATE DUPLICATES THAT ARISE BECAUSE OF EXPANDED LATTICE
          //Check if operation is included in list
          bool contained_in_list = false;
          for (uint i = 0; i < rotation_ops_vec.size(); i++) {
            //Only compare rotations with the same order
            if(aurostd::abs(rotation_ops_vec[i].return_order() - rotation_order) < _ZERO_TOL_) {
              //Use function to compare rotation axes
              if(symmetry_axis_equivalent(L, Linv, point_on_axis, axis_direction, rotation_ops_vec[i].return_point(), rotation_ops_vec[i].return_direction(), tol)) { //DX20190215 - added tol
                if(aurostd::modulus(lattice_vector) < aurostd::modulus(lattice_vectors[i])) {
                  lattice_vectors[i] = lattice_vector;
                }
                contained_in_list = true;
                break;
              }
            }
          }
          //Add to list if not already there
          if(contained_in_list == false) {
            rotation_ops_vec.push_back(candidate_rotation);
            lattice_vectors.push_back(lattice_vector);
            //cerr << tmpchar << endl;
            //cerr << candidate_rotation << endl;
            //cerr << "********************" << endl;
            rotation_ops++;
          }
        }
      }
      cmb = AllCombination42(3, tsize, cmb);
    }
    return rotation_ops_vec;
  }
} //namespace SYM

//[OBSOLETE] // **********************************************************************************************************************
//[OBSOLETE] // closest_point
//[OBSOLETE] // **********************************************************************************************************************
//[OBSOLETE] //Return closest point in cartesian
//[OBSOLETE] xvector<double> closest_point(xvector<double> norm, xmatrix<double> Linv,xmatrix<double> L){
//[OBSOLETE] //Remember that all the symmetry operations of the lattice go through a lattice point. Thus, the representation of any line in the lattice is simple.
//[OBSOLETE] //We do things with reference to "the" origin.
//[OBSOLETE] //Norm is the direction vector of a line (or parametric representation of a line passing through origin) IN CARTESIAN COORDINATES. Linv is the lattice inverse
//[OBSOLETE]   bool LDEBUG=(FALSE || XHOST.DEBUG);
//[OBSOLETE]   xmatrix<double> f2c = trasp(L); //DX
//[OBSOLETE]   xmatrix<double> c2f = inverse(trasp(L)); //DX
//[OBSOLETE]   //DX long long int prec = 1e7;
//[OBSOLETE]   long long int prec = 1e10;
//[OBSOLETE]   //DX double tol = _TOLERANCE_CLOSEST_POINT_;
//[OBSOLETE]   xvector<double> close_point;
//[OBSOLETE]   vector<int> nonzeros;//vector of index where nonzeros are located
//[OBSOLETE]   vector<double> sgns;
//[OBSOLETE]   xvector<double> norm_lattice = f2c*norm;
//[OBSOLETE]   //cerr << "f2c*norm: " << f2c*norm << endl;
//[OBSOLETE]   //xvector<double> norm_lattice = norm*Linv; //(alpha,beta,gamma)
//[OBSOLETE]   int zs=0; //number of zeros
//[OBSOLETE]   //cerr << "norm_lattice: " << norm_lattice << endl;
//[OBSOLETE]   for(int i=1;i<4;i++){
//[OBSOLETE]     if(aurostd::abs(norm_lattice[i])>_SYM_TOL_){
//[OBSOLETE]       nonzeros.push_back(i);
//[OBSOLETE]       sgns.push_back(norm_lattice[i]/aurostd::abs(norm_lattice[i]));
//[OBSOLETE]     }
//[OBSOLETE]     else
//[OBSOLETE]       zs++;
//[OBSOLETE]   }
//[OBSOLETE]   bool cont = true;
//[OBSOLETE]   while(cont == true){
//[OBSOLETE]     if(zs == 3){if(LDEBUG) {cerr << "norm cannot be null: closest_point() in wyckoff_topology_functions" << endl;} xvector<double> null; null[1]=null[2]=null[3]=0.0; return null;}
//[OBSOLETE]     if(zs == 2){
//[OBSOLETE]       close_point(nonzeros[0]) = sgns[0] * 1;
//[OBSOLETE]       //cerr << "2zs" << endl;
//[OBSOLETE]     }
//[OBSOLETE]     if(zs == 1){
//[OBSOLETE]       //cerr << "1zs" << endl;
//[OBSOLETE]       unsigned long long int d1 = cast2int(norm_lattice(nonzeros[0]),prec); //integer*prec
//[OBSOLETE]       unsigned long long int d2 = cast2int(norm_lattice(nonzeros[1]),prec);
//[OBSOLETE]       unsigned long long int GCD = gcd((unsigned long long int) (d1),(unsigned long long int) (d2));
//[OBSOLETE]       //cerr << "d1: " << d1 << endl;
//[OBSOLETE]       //cerr << "d2: " << d2 << endl;
//[OBSOLETE]       //cerr << "GCD: " << GCD << endl;
//[OBSOLETE]       long long int precision_check= GCD % ULLONG_MAX;
//[OBSOLETE]       // === Check if GCD variable has enough memory (if it doesn't precision_check will be negative) === //
//[OBSOLETE]       if(precision_check <0){if(LDEBUG) {cerr << "GCD negative (precision not large enough): closest_point() in wyckoff_topology_functions" << endl;} xvector<double> null; null[1]=null[2]=null[3]=0.0; return null;}
//[OBSOLETE]       close_point(nonzeros[0]) = sgns[0] * double (d1) / double (GCD);
//[OBSOLETE]       close_point(nonzeros[1]) = sgns[1] * double (d2) / double (GCD);
//[OBSOLETE]     }
//[OBSOLETE]     if(zs == 0){
//[OBSOLETE]       //cerr << "0zs" << endl;
//[OBSOLETE]       unsigned long long int d1 = cast2int(norm_lattice(nonzeros[0]),prec);
//[OBSOLETE]       unsigned long long int d2 = cast2int(norm_lattice(nonzeros[1]),prec);
//[OBSOLETE]       unsigned long long int d3 = cast2int(norm_lattice(nonzeros[2]),prec);
//[OBSOLETE]       unsigned long long int GCD = gcd(gcd((unsigned long long int) (d1),(unsigned long long int) (d2)),(unsigned long long int) (d3));
//[OBSOLETE]       long long int precision_check= GCD % ULLONG_MAX;
//[OBSOLETE]       // === Check if GCD variable has enough memory (if it doesn't precision_check will be negative) === //
//[OBSOLETE]       if(precision_check <0){if(LDEBUG) {cerr << "GCD negative (precision not large enough): closest_point() in wyckoff_topology_functions" << endl;} xvector<double> null; null[1]=null[2]=null[3]=0.0; return null;}
//[OBSOLETE]       close_point(nonzeros[0]) = sgns[0] * double (d1) / double (GCD);
//[OBSOLETE]       close_point(nonzeros[1]) = sgns[1] * double (d2) / double (GCD);
//[OBSOLETE]       close_point(nonzeros[2]) = sgns[2] * double (d3) / double (GCD);
//[OBSOLETE]     }
//[OBSOLETE]     //Now have closest point in terms of lattice.
//[OBSOLETE]    double largest_val =  infnorm<double>(close_point);
//[OBSOLETE]     if(largest_val/prec > _SYM_TOL_){
//[OBSOLETE]       prec *= 10;
//[OBSOLETE]     }
//[OBSOLETE]     else {
//[OBSOLETE]       cont = false;
//[OBSOLETE]     }
//[OBSOLETE]   }
//[OBSOLETE]   //Transform into cartesian coordinates for output:
//[OBSOLETE]   close_point = close_point*L;
//[OBSOLETE]   for(int i=1;i<4;i++){
//[OBSOLETE]     if(aurostd::abs(close_point[i]) < 1/double(prec)){
//[OBSOLETE]       close_point[i] = 0.0;
//[OBSOLETE]     }
//[OBSOLETE]   }
//[OBSOLETE]   return close_point;
//[OBSOLETE] }

//Number Theory functions:
//[OBSOLETE]// **********************************************************************************************************************
//[OBSOLETE]// cast2int
//[OBSOLETE]// **********************************************************************************************************************
//[OBSOLETE]long long int cast2int(double d, long long int prec) {
//[OBSOLETE]  //double tol = 1e-10;
//[OBSOLETE]  double tol = _SYM_TOL_;
//[OBSOLETE]  long long int a = 0;
//[OBSOLETE]  if(aurostd::abs(ceil(aurostd::abs(d)) - aurostd::abs(d)) < tol) {
//[OBSOLETE]    a = (long long int)ceil(aurostd::abs(d)) * prec;
//[OBSOLETE]  } else if(aurostd::abs(aurostd::abs(d) - floor(aurostd::abs(d))) < tol) {
//[OBSOLETE]    a = (long long int)floor(aurostd::abs(d)) * prec;
//[OBSOLETE]  } else {
//[OBSOLETE]    a = (long long int)(aurostd::abs(d) * prec);
//[OBSOLETE]  }
//[OBSOLETE]  return a;
//[OBSOLETE]}

// **********************************************************************************************************************
// gcdD (Dijkstra's Algorithm)
// **********************************************************************************************************************
//Dijkstra's Algorithm
//DX20191202 [OBSOLETE]namespace SYM {
//DX20191202 [OBSOLETE]  int gcdD(int m, int n) {
//DX20191202 [OBSOLETE]    if(m == 1 || n == 1)
//DX20191202 [OBSOLETE]      return 1;
//DX20191202 [OBSOLETE]    if(m == n)
//DX20191202 [OBSOLETE]      return m;
//DX20191202 [OBSOLETE]    else if(m > n)
//DX20191202 [OBSOLETE]      return gcdD(m - n, n);
//DX20191202 [OBSOLETE]    else
//DX20191202 [OBSOLETE]      return gcdD(m, n - m);
//DX20191202 [OBSOLETE]  }
//DX20191202 [OBSOLETE]} //namespace SYM

//ASSUME POSITIVE INPUTS

// **********************************************************************************************************************
// gcd
// **********************************************************************************************************************
//DX20191202 [OBSOLETE]namespace SYM {
//DX20191202 [OBSOLETE]  unsigned long long int gcd(unsigned long long int m, unsigned long long int n) {
//DX20191202 [OBSOLETE]    //cerr << "GCD algorithm " << m << " " << n <<  endl;
//DX20191202 [OBSOLETE]    //double tol = 1e-6;
//DX20191202 [OBSOLETE]    double tol = _ZERO_TOL_;
//DX20191202 [OBSOLETE]    //make largest input 'm'
//DX20191202 [OBSOLETE]    if(n > m) {
//DX20191202 [OBSOLETE]      unsigned long long int tmp = n;
//DX20191202 [OBSOLETE]      n = m;
//DX20191202 [OBSOLETE]      m = tmp;
//DX20191202 [OBSOLETE]    }
//DX20191202 [OBSOLETE]    //cast to INT
//DX20191202 [OBSOLETE]    unsigned long long int a1 = (unsigned long long int)m;
//DX20191202 [OBSOLETE]    unsigned long long int a2 = (unsigned long long int)n;
//DX20191202 [OBSOLETE]    unsigned long long int b;
//DX20191202 [OBSOLETE]
//DX20191202 [OBSOLETE]    //Check if one of the inputs is unity:
//DX20191202 [OBSOLETE]    if(aurostd::abs(m - 1) < tol || aurostd::abs(n - 1) < tol) {
//DX20191202 [OBSOLETE]      return 1.0;
//DX20191202 [OBSOLETE]    }
//DX20191202 [OBSOLETE]    //Check if the inputs are equal:
//DX20191202 [OBSOLETE]    long double dividend = aurostd::abs(m / n);
//DX20191202 [OBSOLETE]    if(aurostd::abs(dividend - floor(dividend)) < tol || aurostd::abs(dividend - ceil(dividend)) < tol)
//DX20191202 [OBSOLETE]      return n;
//DX20191202 [OBSOLETE]    //Procede with algorithm:
//DX20191202 [OBSOLETE]    /////////////////////
//DX20191202 [OBSOLETE]    //NON RECURSIVE FUNCTION:
//DX20191202 [OBSOLETE]    else {
//DX20191202 [OBSOLETE]      while ((b = a1 % a2) != 0) {
//DX20191202 [OBSOLETE]	a1 = a2;
//DX20191202 [OBSOLETE]	a2 = b;
//DX20191202 [OBSOLETE]      }
//DX20191202 [OBSOLETE]    }
//DX20191202 [OBSOLETE]    return (unsigned long long int)a2;
//DX20191202 [OBSOLETE]    ////////////////////
//DX20191202 [OBSOLETE]
//DX20191202 [OBSOLETE]    //else if(aurostd::abs(m/n-1) > tol)
//DX20191202 [OBSOLETE]    //  return gcd(m-n, n);
//DX20191202 [OBSOLETE]    //else
//DX20191202 [OBSOLETE]    //  return gcd(m, n-m);
//DX20191202 [OBSOLETE]  }
//DX20191202 [OBSOLETE]} //namespace SYM
//DX20191202 [OBSOLETE]
//DX20191202 [OBSOLETE]// **********************************************************************************************************************
//DX20191202 [OBSOLETE]// gcd (Euclid's Algorithm)
//DX20191202 [OBSOLETE]// **********************************************************************************************************************
//DX20191202 [OBSOLETE]//Euclid's algorithm
//DX20191202 [OBSOLETE]namespace SYM {
//DX20191202 [OBSOLETE]  long long int gcd(long long int x, long long int y) {
//DX20191202 [OBSOLETE]    if(y > x) {
//DX20191202 [OBSOLETE]      int tmp = x;
//DX20191202 [OBSOLETE]      x = y;
//DX20191202 [OBSOLETE]      y = tmp;
//DX20191202 [OBSOLETE]    }
//DX20191202 [OBSOLETE]    //  cerr << x << " " << y << endl;
//DX20191202 [OBSOLETE]    if(y == 0)
//DX20191202 [OBSOLETE]      return x;
//DX20191202 [OBSOLETE]    else if(x >= y && y > 0)
//DX20191202 [OBSOLETE]      return gcd(y, x % y);
//DX20191202 [OBSOLETE]    else
//DX20191202 [OBSOLETE]      return -1;
//DX20191202 [OBSOLETE]  }
//DX20191202 [OBSOLETE]} //namespace SYM

// **********************************************************************************************************************
// //DX20210422 [OBSOELTE] allsame -> aurostd::identical_entries() //DX20210422
// **********************************************************************************************************************

// **********************************************************************************************************************
// ReturnITCGenShift
// **********************************************************************************************************************
//Returns a vector of pure shifts of generator operations from the ITC tables for a given spacegroup.
namespace SYM {
  vector<xvector<double> > ReturnITCGenShift(int spacegroupnum, string axis_cell) {
    SymmetryInformationITC ITC_sym_info; //DX20190215
    ITC_sym_info.initsgs(axis_cell); //DX20190215
    //DX20190215 [OBSOLETE]#ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
    //DX20190215 [OBSOLETE]    extern thread_local vector<string> gl_sgs;
    //DX20190215 [OBSOLETE]#else
    //DX20190215 [OBSOLETE]    extern vector<string> gl_sgs;
    //DX20190215 [OBSOLETE]#endif
    spacegroupnum = spacegroupnum - 1;
    vector<int> genlocations = generatorindex(spacegroupnum + 1);
    vector<xvector<double> > out;
    int s = 0;
    vector<int> wint = get_multiplicities(ITC_sym_info.gl_sgs[spacegroupnum]); //DX20190215
    vector<vector<vector<string> > > tmpvvvstring = get_wyckoff_pos(ITC_sym_info.gl_sgs[spacegroupnum], wint[1], false);  //false = do not get centered points
    vector<vector<vector<vector<sdouble> > > > tmpvvvsd = convert_wyckoff_pos_sd(tmpvvvstring);

    //cerr << gl_sgs[spacegroupnum] << endl;
    //cerr << "WYCKOFF POSTIONS: " << endl;
    //print_wyckoff_pos(tmpvvvstring);
    //throw aurostd::xerror(__AFLOW_FILE__,XPID+"SYM::ReturnITCGenShift():","Throw for debugging purposes.",_GENERIC_ERROR_);

    //for(int k=0;k<tmpvvvsd[s].size();k++){
    //  for(int j=0;j<tmpvvvsd[s][k].size();j++){
    //    xb();
    //    for(int i=0;i<tmpvvvsd[s][k][j].size();i++){
    //	cerr << tmpvvvsd[s][k][j][i].dbl;
    //	cerr << tmpvvvsd[s][k][j][i].chr << " ";
    //    }
    //  }
    //}

    for (uint k = 0; k < tmpvvvsd[s].size(); k++) {
      //DX20210421 [OBSOLETE] if(intinvec(genlocations, k + 1))
      if(aurostd::WithinList(genlocations, (int)k + 1)) { //DX20210421
        xvector<double> oneshift;
        for (uint j = 0; j < tmpvvvsd[s][k].size(); j++) {
          double filler = 0;
          for (uint i = 1; i < tmpvvvsd[s][k][j].size(); i++) {
            //cerr <<"dbl: " <<  tmpvvvsd[s][k][j][i].dbl << endl;
            //cerr << "chr: " <<  tmpvvvsd[s][k][j][i].chr << endl;
            char tmpchr = tmpvvvsd[s][k][j][i].chr;
            if(tmpchr == '\0') {
              double tmpdbl = tmpvvvsd[s][k][j][i].dbl;
              filler = filler + tmpdbl;
            }
          }
          oneshift(j + 1) = filler;
        }
        out.push_back(oneshift);
        oneshift.clear();
      }
    }
    return out;
  }
} //namespace SYM

// **********************************************************************************************************************
// ReverseSpaceGroup
// **********************************************************************************************************************
namespace SYM {
  string ReverseSpaceGroup(vector<string> num) {
    stringstream oss;
    int spacegroupnum = atoi(num[2].c_str()) - 1;
    int m = atoi(num[3].c_str());  // multiplicity of wyckoff site (eg, 48 j 1 (TWICE))
    int n = 1;
    int l = 1;

    if(num.size() == 6) {
      n = atoi(num[4].c_str());  // which multiplicity group
      l = atoi(num[5].c_str());  // which sub group within multiplicity group (first, second, third, etc)
    }
    if(l < 1) {
      throw aurostd::xerror(__AFLOW_FILE__,XPID+"SYM::ReverseSpaceGroup():","Sub multiplicity groups start at 1.",_INPUT_ILLEGAL_);
    }

    string axis_cell = "";
    SymmetryInformationITC ITC_sym_info; //DX20190215
    ITC_sym_info.initsgs(axis_cell); //DX20190215
    //DX20190215 [OBSOLETE] #ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
    //DX20190215 [OBSOLETE]     extern thread_local vector<string> gl_sgs;
    //DX20190215 [OBSOLETE] #else
    //DX20190215 [OBSOLETE]     extern vector<string> gl_sgs;
    //DX20190215 [OBSOLETE] #endif
    xmatrix<double> L = spacegroup_to_lattice(spacegroupnum + 1);

    get_random_double(.1, .9);
    double r1;
    double r2;
    double r3;
    vector<int> wint_all = get_multiplicities(ITC_sym_info.gl_sgs[spacegroupnum]); //DX20190215
    vector<int> wint; //DX20210422
    for(uint i=0;i<wint_all.size();i++){ //DX20210422
      if(!aurostd::WithinList(wint,wint_all[i])){ wint.push_back(wint_all[i]); } //DX20210422
    } //DX20210422
    //TITLE
    oss << "REVERSE SPACE GROUP: " << spacegroupnum + 1 << " using WYCKOFF multiplicity:" << wint[n] << "-" << l << " " << endl;
    oss << "1" << endl;
    oss << L << endl;
    vector<vector<vector<string> > > tmpvvvstring = get_wyckoff_pos(ITC_sym_info.gl_sgs[spacegroupnum], wint[n]);
    //digest_wyckoff_pos(gl_sgs[spacegroupnum]);
    vector<vector<vector<vector<sdouble> > > > tmpvvvsd = convert_wyckoff_pos_sd(tmpvvvstring);
    int s = l - 1;

    //cerr << gl_sgs[spacegroupnum] << endl;
    //print_wyckoff_pos(tmpvvvstring);

    vector<_atom> atoms;
    int atomlabel = -1;
    int countold = 0;
    int total = 0;
    for (int d = 0; d < m; d++) {
      atomlabel++;
      stringstream ss;
      ss << atomlabel;
      string typestring = ss.str();
      //r1 = get_random_double(.1,.4);
      //r2 = get_random_double(.5,.6);
      //r3 = get_random_double(.7,.8);
      r1 = 1 / ((double)d + 1);
      r2 = 1 / ((double)d + 1.5);
      r3 = 1 / ((double)d + 1.8);
      //r1 = .4;
      //r2 = .2;
      //r3 = .3;
      for (uint k = 0; k < tmpvvvsd[s].size(); k++) {
        _atom tmp;
        tmp.name = typestring;
        for (uint j = 0; j < tmpvvvsd[s][k].size(); j++) {
          //xb();
          for (uint i = 0; i < tmpvvvsd[s][k][j].size(); i++) {
            //cerr <<"dbl: " <<  tmpvvvsd[s][k][j][i].dbl << endl;
            //cerr << "chr: " <<  tmpvvvsd[s][k][j][i].chr << endl;
            char tmpchr = tmpvvvsd[s][k][j][i].chr;
            double tmpdbl = tmpvvvsd[s][k][j][i].dbl;

            if(tmpchr == 'x') {
              tmp.fpos(j + 1) += tmpdbl * r1;
            } else if(tmpchr == 'y') {
              tmp.fpos(j + 1) += tmpdbl * r2;
            } else if(tmpchr == 'z') {
              tmp.fpos(j + 1) += tmpdbl * r3;
            } else {
              tmp.fpos(j + 1) += tmpdbl;
            }
          }
        }
        atoms.push_back(tmp);
      }
      total += atoms.size();
      oss << atoms.size() - countold << " ";
      countold = atoms.size();
    }
    oss << endl;
    oss << "Direct(" << total << ")" << endl;
    ostringstream css;
    for (uint i = 0; i < atoms.size(); i++) {
      oss << fixed << setprecision(15) << atoms[i].fpos(1) << " " << atoms[i].fpos(2) << " " << atoms[i].fpos(3) << " " << atoms[i].name << endl;
    }
    return oss.str();
  }
} //namespace SYM

// **********************************************************************************************************************
// points_equivalent
// **********************************************************************************************************************
//L_inv is the inverse of the lattice under question (makes loops quicker not doing inverse n times).
namespace SYM {
  bool points_equivalent(xmatrix<double>& c2f, xmatrix<double>& f2c, xvector<double> P1, xvector<double> P2, xvector<double>& lattice_vector, double& radius, bool& skew, double& tol) { //DX20190215
    xvector<double> P1_fpos = c2f * P1;
    xvector<double> P2_fpos = c2f * P2;

    xvector<double> fdiff = P1_fpos - P2_fpos;
    double orig_cdiff = aurostd::modulus(f2c * fdiff);
    xvector<int> ijk;
    //DX20190613 [OBSOLETE] bool neighbor_restriction = false;
    double min_dist = 1e9; //DX20190613

    if(skew) {
      xmatrix<double> lattice = trasp(f2c); //DX20190613
      xvector<double> min_cdiff = SYM::minimizeDistanceCartesianMethod(P1, P2, lattice, ijk); //DX20190613
      fdiff = c2f * min_cdiff; //DX20190613
      min_dist = aurostd::modulus(min_cdiff); //DX20190613
      //DX20190613 [OBSOLETE] SYM::minimizeCartesianDistance(P1, P2, fdiff, c2f, f2c, ijk, neighbor_restriction, tol);
    } else {
      fdiff = SYM::minimizeDistanceFractionalMethod(fdiff,ijk); //DX20190613
      min_dist = aurostd::modulus(f2c*fdiff); //DX20190613
      //DX20190613 [OBSOLETE] SYM::PBC(fdiff, ijk, neighbor_restriction);
    }
    for (uint i = 1; i < 4; i++) {
      if(aurostd::abs(ijk(i)) > 3 && (orig_cdiff - radius) > _ZERO_TOL_) {
        return false;
      }
    }
    xvector<double> tmp_fdiff;
    for (uint i = 1; i < 4; i++) {
      tmp_fdiff(i) = fdiff(i) + (double)ijk(i);
    }
    lattice_vector = f2c * tmp_fdiff;
    return (min_dist < tol); //DX20190613
    //DX20190613 [OBSOLETE] return (aurostd::modulus(f2c * fdiff) < tol);
  }
} //namespace SYM

// **********************************************************************************************************************
// is_lattice_point
// **********************************************************************************************************************
namespace SYM {
  bool is_lattice_point(xmatrix<double> L, xvector<double> point, xvector<double>& lattice_vector, double& radius, bool& skew, double& tol) { //DX20190215
    xmatrix<double> f2c = trasp(L);
    xmatrix<double> c2f = inverse(trasp(L));
    xvector<double> tmp = c2f * point;

    xvector<double> origin;
    origin(1) = 0.0;
    origin(2) = 0.0;
    origin(3) = 0.0;
    xvector<double> fdiff = origin - tmp;
    //xvector<double> fdiff = tmp - origin;
    double orig_cdiff = aurostd::modulus(point);
    xvector<int> ijk;
    double min_dist = 1e9; //DX20190613
    //DX20190613 [OBSOLETE] bool neighbor_restriction = false;
    if(skew) {
      xvector<double> min_cdiff = SYM::minimizeDistanceCartesianMethod(origin, point, L, ijk); //DX20190613
      fdiff = c2f * min_cdiff; //DX20190613
      min_dist = aurostd::modulus(min_cdiff); //DX20190613
      //DX20190613 [OBSOLETE] SYM::minimizeCartesianDistance(origin, point, fdiff, c2f, f2c, ijk, neighbor_restriction, tol); //DX20190215
    } else {
      fdiff = SYM::minimizeDistanceFractionalMethod(fdiff, ijk); //DX20190613
      min_dist = aurostd::modulus(f2c * fdiff); //DX20190613
      //DX20190613 [OBSOLETE] SYM::PBC(fdiff, ijk, neighbor_restriction);
    }
    xvector<double> tmp_fdiff;
    for (uint i = 1; i < 4; i++) {
      if(aurostd::abs(ijk(i)) > 3 && (orig_cdiff - radius) > _ZERO_TOL_) {
        return false;
      }
      tmp_fdiff(i) = fdiff(i) + (double)ijk(i);
    }
    lattice_vector = f2c * tmp_fdiff;
    return (min_dist < tol); //DX20190215
    //DX20190613 [OBSOLETE] return (aurostd::modulus(f2c * fdiff) < tol); //DX20190215
  }
} //namespace SYM

// **********************************************************************************************************************
// screw_equivalent
// **********************************************************************************************************************
//Check if two screws are the same (assuming they both exist in the same unit cell)
namespace SYM {
  bool screw_equivalent(Screw S1, Screw S2) {
    bool same = false;
    double tol = 1e-6;
    xvector<double> dir1 = S1.return_direction();
    xvector<double> dir2 = S2.return_direction();
    dir1/=aurostd::modulus(dir1); //DX20190905 - SYM::normalize() -> divide by aurostd::modulus()
    dir2/=aurostd::modulus(dir2); //DX20190905 - SYM::normalize() -> divide by aurostd::modulus()
    if(dir1 == dir2) {
      if(aurostd::modulus(aurostd::vector_product(S1.return_point() - S2.return_point(), S1.return_direction())) < tol) { //DX20190905 - SYM::CrossPro() -> aurostd::vector_product()
        same = true;
      }
    }
    return same;
  }
} //namespace SYM

// **********************************************************************************************************************
// symmetry_axis_equivalent
// **********************************************************************************************************************
namespace SYM {
  bool symmetry_axis_equivalent(xmatrix<double> L, xmatrix<double> Linv, xvector<double> P1, xvector<double> N1, xvector<double> P2, xvector<double> N2, double& tol) { //DX20190215
    //DX20190215 [OBSOLETE] double tol = _SYM_TOL_;
    bool same = true;
    int count = 0;
    //First the axis directions must be parallel
    //Normalize the axis/plane normal directions
    N1 = (1 / aurostd::modulus(N1)) * N1;
    N2 = (1 / aurostd::modulus(N2)) * N2;
    if(aurostd::modulus((N1 - N2) * L) < tol || aurostd::modulus((-N1 - N2) * L) < tol) {
      xvector<double> T;
      xvector<double> tmp;
      //double d1,d2,d3;
      int check;
      vector<double> nonzeros;

      //(p2+a*N)-p1 = (ijk)*L
      //if a=0 this can be solved using L_inv: (p2-p1)*L_inv must be integers
      //otherwise have to prove that an 'a' exists where ((p2+a*N)-p1)*L_inv = integers
      tmp = (P2 - P1) * Linv;
      for (int i = 1; i < 4; i++) {
        if(tmp(i) - floor(tmp(i)) < tol || ceil(tmp(i)) - tmp(i) < tol)
          count++;
      }
      if(count == 3) {
        return same;
      } else {
        for (int i = -2; i <= 2; i++) {
          for (int j = -2; j <= 2; j++) {
            for (int k = -2; k <= 2; k++) {
              nonzeros.clear();
              check = 0;
              T(1) = i;
              T(2) = j;
              T(3) = k;
              //cerr << T << endl;
              tmp = T * L + P1 - P2;
              for (int ii = 1; ii < 4; ii++) {
                if(aurostd::abs(tmp(ii)) < tol && aurostd::abs(N1(ii)) < tol)
                  check++;
                if(aurostd::abs(N1(ii)) > tol) {
                  nonzeros.push_back(tmp(ii) / N1(ii));
                }
              }
              if(aurostd::identical(nonzeros,_ZERO_TOL_) && check + nonzeros.size() == 3) { //DX20210422 - SYM::allsame() -> aurostd::identical_entries()
                same = true;
                return same;
              } else {
                same = false;
                //	      return same;
              }
            }
          }
        }
      }
      return same;
    } else {
      same = false;
      return same;
    }
  }
} //namespace SYM

// **********************************************************************************************************************
// random_point
// **********************************************************************************************************************
namespace SYM {
  xvector<double> random_point() {
    xvector<double> out;
    srand(time(0));
    out(1) = double(rand() % 50 + 100) / 1000;
    out(2) = double(rand() % 50 + 151) / 1000;
    out(3) = double(rand() % 50 + 202) / 1000;
    return out;
  }
} //namespace SYM

// **********************************************************************************************************************
// add_3d_point
// **********************************************************************************************************************
namespace SYM {
  void add_3d_point(vector<xvector<double> >& points, double x, double y, double z) {
    xvector<double> tmp;
    tmp(1) = x;
    tmp(2) = y;
    tmp(3) = z;
    points.push_back(tmp);
  }
} //namespace SYM

// **********************************************************************************************************************
// get_mod_angle_tensor
// **********************************************************************************************************************
namespace SYM {
  xmatrix<double> get_mod_angle_tensor(vector<xvector<double> >& points) {
    //RETURNS A MATRIX WITH VECTOR MAGNITUDES ALONG DIAGONAL AND ANGLES BETWEEN VECTORS ON OFF DIAGONAL PARTS
    //Check that all points are same dimension
    int size = points[0].urows;
    for (uint i = 0; i < points.size(); i++) {
      if(points[i].urows != size) {
        throw aurostd::xerror(__AFLOW_FILE__,XPID+"SYM::get_mod_angle_tensor():","Points must have the same dimension.",_VALUE_RANGE_);
      }
    }
    //Choose first point as origin
    vector<xvector<double> > points_referenced_to_first;
    for (uint i = 1; i < points.size(); i++) {
      points_referenced_to_first.push_back(points[i] - points[0]);
    }
    size = points_referenced_to_first.size();

    xmatrix<double> metric_tensor(size, size);

    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        metric_tensor(i + 1, j + 1) = aurostd::scalar_product(points_referenced_to_first[i], points_referenced_to_first[j]);
        if(i != j) {
          metric_tensor(i + 1, j + 1) = acos(metric_tensor(i + 1, j + 1) / (aurostd::modulus(points_referenced_to_first[i]) * aurostd::modulus(points_referenced_to_first[j])));
        }
        if(i == j) {
          metric_tensor(i + 1, j + 1) = sqrt(metric_tensor(i + 1, j + 1));
        }
      }
    }
    return metric_tensor;
  }
} //namespace SYM

// **********************************************************************************************************************
// discern_rot_sym
// **********************************************************************************************************************
//returns char for potential rotation types: 'a' = 3fold 'b' = 4fold 'c' = 6fold
namespace SYM {
  char discern_rot_sym(xmatrix<double> m, double& tol) {  //m is metric_tensor //DX20190215
    double two_fold = Pi_r / 2.0;
    double four_fold = Pi_r / 4.0;
    double six_fold_120 = (2.0 * Pi_r) / (3.0);
    double six_fold_60 = (Pi_r) / (6.0);
    char out = 'z';
    //Check 3-fold
    if(aurostd::abs(m(1, 1) - m(2, 2)) < tol && ((m(1, 1) + m(2, 2)) / 2) * sin(aurostd::abs(m(1, 2) - Pi_r / 3)) < tol) { //DX20190215
      out = 'a';
      return out;
    }

    //Check 4-fold
    if(aurostd::abs(m(1, 1) - m(2, 2)) < tol) {
      if(SYM::checkAngle(m(1, 1), m(2, 2), m(1, 2), two_fold, tol)) { //DX20190215
        out = 'b';
        return out;
      }
    }
    if(((m(1, 1) / m(2, 2)) - sqrt(2) < tol || (m(1, 1) / m(2, 2)) - (1 / sqrt(2)) < tol) && SYM::checkAngle(m(1, 1), m(2, 2), m(1, 2), four_fold, tol)) { //DX20190215
      out = 'b';
      return out;
    }
    //Check 6-fold
    if(aurostd::abs(m(1, 1) - m(2, 2)) < tol) {
      if(SYM::checkAngle(m(1, 1), m(2, 2), m(1, 2), six_fold_120, tol)) { //DX20190215
        out = 'c';
        return out;
      }
    }
    if(((m(1, 1) / m(2, 2)) - sqrt(3) < tol || (m(1, 1) / m(2, 2)) - (1 / sqrt(3)) < tol) && SYM::checkAngle(m(1, 1), m(2, 2), m(1, 2), six_fold_60, tol)) { //DX20190215
      out = 'c';
      return out;
    }

    else
      return out;
  }
} //namespace SYM

#endif
