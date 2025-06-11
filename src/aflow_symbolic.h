// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// Written by David Hicks (DX) - 2020

#ifndef _AFLOW_SYMBOLIC_H_
#define _AFLOW_SYMBOLIC_H_

#if USE_SYMBOLIC_SOURCE //DX20200831 - defined in aflow.h

#include "aflow.h"
#include "aflow_symmetry_spacegroup.h"
#include "aflow_compare_structure.h"
#include "SYMBOLICCPLUSPLUS/symbolic_main.h"

//// toggle symbolic math
//#define COMPILE_SYMBOLIC


// Symbolic math variables
#define _SYMBOLIC_ZERO_ symbolic::Symbolic(0)
#define _ANRL_LATTICE_VARIABLES_  "a,b,c,alpha,beta,gamma,cx,cy,cz"
#define _ANRL_TRIG_VARIABLES_  "sin,cos,tan,sec,csc,cot"
#define _ANRL_WYCKOFF_VARIABLES_  "x,y,z"

//SYMBOLIC

template <typename T> char const* get_type(T const& object);

// ---------------------------------------------------------------------------
// mirrors wyckoffsite_ITC, requires SymbolicC++ 
// instead of adding Symbolic equations to wyckoffsite_ITC (don't want to put 
// large dependency in aflow.h)
struct SymbolicWyckoffSite {
  xvector<double> coord;
  uint index;
  string type;
  string wyckoffSymbol;
  string letter;
  string site_symmetry;
  uint multiplicity;
  double site_occupation;
  vector<symbolic::Symbolic> equations;
  uint parameter_index;
};

SymbolicWyckoffSite initializeSymbolicWyckoffSite(const wyckoffsite_ITC& Wyckoff);
void substituteVariableWithParameterDesignation(vector<SymbolicWyckoffSite>& Wyckoff_sites_symbolic);
void substituteVariableWithParameterDesignation(SymbolicWyckoffSite& Wyckoff_symbolic);

// ---------------------------------------------------------------------------
// "pure" symbolic functions 
namespace symbolic {
  Symbolic string2symbolic(const string& str);
  bool isEqual(const Symbolic& a, const Symbolic& b);
  bool isEqualVector(const Symbolic& a_vec, const Symbolic& b_vec);
  vector<vector<string> > matrix2VectorVectorString(const Symbolic& lattice);
  Symbolic BringInCell(const Symbolic& vec_in, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0);
}

// ---------------------------------------------------------------------------
// ANRL symbolic functions 
namespace anrl {
  symbolic::Symbolic SymbolicANRLPrimitiveLattices(const string& lattice_and_centering, const char& space_group_letter);
  vector<symbolic::Symbolic> equations2SymbolicEquations(const vector<vector<string> >& equations);
  symbolic::Symbolic cartesian2lattice(const symbolic::Symbolic& lattice, const symbolic::Symbolic& cartesian_coordinate);
  symbolic::Symbolic getXYZ2LatticeTransformation(const string& lattice_and_centering);
  vector<symbolic::Symbolic> getEquationsForCenteredLattices(const string& lattice_and_centering,
      const symbolic::Symbolic& lattice,
      const vector<symbolic::Symbolic>& conventional_equations);
  vector<symbolic::Symbolic> convertEquations2FractionalEquations(const string& lattice_and_centering,
      const symbolic::Symbolic& lattice,
      const vector<symbolic::Symbolic> conventional_equations);
  void addSymbolicEquation2Atoms(const vector<symbolic::Symbolic>& equations, deque<_atom>& atoms, bool isfpos=true);
}
#endif

#endif
