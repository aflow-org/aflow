// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// Written by David Hicks (DX) - 2020

#ifndef _AFLOW_SYMBOLIC_H_
#define _AFLOW_SYMBOLIC_H_

#include <deque>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_defs.h"
#include "extern/SYMBOLICCPLUSPLUS/symbolic.h"
#include "structure/aflow_xatom.h"

#if USE_SYMBOLIC_SOURCE // DX20200831 - defined in aflow.h
#include "extern/SYMBOLICCPLUSPLUS/symbolic_main.h" // NOLINT / Needed because the way symbolic is included is a mess

//// toggle symbolic math
// #define COMPILE_SYMBOLIC

// Symbolic math variables
#define _SYMBOLIC_ZERO_ symbolic::Symbolic(0)
#define _ANRL_LATTICE_VARIABLES_ "a,b,c,alpha,beta,gamma,cx,cy,cz"
#define _ANRL_TRIG_VARIABLES_ "sin,cos,tan,sec,csc,cot"
#define _ANRL_WYCKOFF_VARIABLES_ "x,y,z"

// SYMBOLIC

template <typename T> const char* get_type(const T& object);

// ---------------------------------------------------------------------------
// mirrors wyckoffsite_ITC, requires SymbolicC++
// instead of adding Symbolic equations to wyckoffsite_ITC (don't want to put
// large dependency in aflow.h)
struct SymbolicWyckoffSite {
  aurostd::xvector<double> coord;
  uint index;
  std::string type;
  std::string wyckoffSymbol;
  std::string letter;
  std::string site_symmetry;
  uint multiplicity;
  double site_occupation;
  std::vector<symbolic::Symbolic> equations;
  uint parameter_index;
};

SymbolicWyckoffSite initializeSymbolicWyckoffSite(const wyckoffsite_ITC& Wyckoff);
void substituteVariableWithParameterDesignation(std::vector<SymbolicWyckoffSite>& Wyckoff_sites_symbolic);
void substituteVariableWithParameterDesignation(SymbolicWyckoffSite& Wyckoff_symbolic);

// ---------------------------------------------------------------------------
// "pure" symbolic functions
namespace symbolic {
  Symbolic string2symbolic(const std::string& str);
  bool isEqual(const Symbolic& a, const Symbolic& b);
  bool isEqualVector(const Symbolic& a_vec, const Symbolic& b_vec);
  std::vector<std::vector<std::string>> matrix2VectorVectorString(const Symbolic& lattice);
  Symbolic BringInCell(const Symbolic& vec_in, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0);
} // namespace symbolic

// ---------------------------------------------------------------------------
// ANRL symbolic functions
namespace anrl {
  symbolic::Symbolic SymbolicANRLPrimitiveLattices(const std::string& lattice_and_centering, const char& space_group_letter);
  std::vector<symbolic::Symbolic> equations2SymbolicEquations(const std::vector<std::vector<std::string>>& equations);
  symbolic::Symbolic cartesian2lattice(const symbolic::Symbolic& lattice, const symbolic::Symbolic& cartesian_coordinate);
  symbolic::Symbolic getXYZ2LatticeTransformation(const std::string& lattice_and_centering);
  std::vector<symbolic::Symbolic> getEquationsForCenteredLattices(const std::string& lattice_and_centering, const symbolic::Symbolic& lattice, const std::vector<symbolic::Symbolic>& conventional_equations);
  std::vector<symbolic::Symbolic> convertEquations2FractionalEquations(const std::string& lattice_and_centering, const symbolic::Symbolic& lattice, const std::vector<symbolic::Symbolic> conventional_equations);
  void addSymbolicEquation2Atoms(const std::vector<symbolic::Symbolic>& equations, std::deque<_atom>& atoms, bool isfpos = true);
} // namespace anrl
#endif

#endif
