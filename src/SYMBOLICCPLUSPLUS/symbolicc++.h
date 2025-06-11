//   SymbolicC++ : An object oriented computer algebra system written in C++
//
//   Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb
//
//   This library is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This library is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License along
//   with this program; if not, write to the Free Software Foundation, Inc.,
//   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

// symbolicc++.h

#include "../AUROSTD/aurostd.h" //DX20200817 - to include std namespace libraries

// normal include headers
#include <iostream>
#include <iterator>
#include <list>
#include "cloning.h"
#include "identity.h"

// phased include headers
// according to class hierarchy
// DX20200612 - removed class hierarchy (one-level, i.e., removed symbolic directory and moved everything to the same level; see README_AFLOW_SYMBOLIC.TXT)

// DX: Note, the order of includes is important (dependencies and forward declarations), do not change the order

#include "symbolic.h"  // SymbolicInterface, Symbolic ...
#include "equation.h"  //   Equation : CloningSymbolicInterface
#include "number.h"    //   Numeric  : CloningSymbolicInterface
#include "matrix.h"
#include "product.h"   //   Product  : CloningSymbolicInterface
#include "sum.h"       //   Sum      : CloningSymbolicInterface
#include "symbol.h"    //   Symbol   : CloningSymbolicInterface
#include "functions.h" //     Sin    : Symbol ...
#include "symmatrix.h" //   SymbolicMatrix : CloningSymbolicInterface
#include "symerror.h"  //   SymbolicError  : CloningSymbolicInterface
#include "constants.h" 
#include "integrate.h"
#include "solve.h"

//Extra headers (not included in original source code, i.e., these functions were not built by default) - DX20200824
#include "array.h"
#include "derive.h"
#include "matnorm.h"
#include "multinomial.h"
#include "polynomial.h"
#include "quatern.h"
#include "rational.h"
#include "vecnorm.h"
#include "vector.h"
#include "verylong.h"

#ifndef SYMBOLIC_CPLUSPLUS
#define SYMBOLIC_CPLUSPLUS

// Include the relevant classes in 3 phases
//  Phase 1 ensures that every class has a forward
//  declaration for use in phase 2.
//  Phase 2 ensures that every constructor and method
//  has a forward declaration for use in phase 3.
// 1. Forward declarations: class X;
// 2. Declaraions:          class X { ... };
// 3. Definitions:          X::X() ...

// This overcomes mutual recursion in dependencies,
// for example class Sum needs class Product
// and class Product needs class Sum.

// forward declarations of all classes first
//DX20200831 [OBSOLETE] #define SYMBOLIC_FORWARD
//DX20200831 [OBSOLETE] #include "symbolicc++.h"
//DX20200831 [OBSOLETE] #undef  SYMBOLIC_FORWARD

namespace symbolic{ //DX20200625
  typedef std::list<Equation> Equations;
  typedef std::list<Equations> PatternMatches;
} //namespace symbolic //DX20200625

// declarations of classes without definitions
#define SYMBOLIC_DECLARE
#include "symbolicc++.h"
#undef  SYMBOLIC_DECLARE

// declarations for non-member functions
// also used in definition phase for clarity

namespace symbolic{ //DX20200625
  Symbolic expand(const SymbolicInterface&);
  ostream &operator<<(ostream &,const Symbolic &);
  ostream &operator<<(ostream &,const Equation &);
  Symbolic operator+(const Symbolic &);
  Symbolic operator+(const Symbolic &,const Symbolic &);
  Symbolic operator+(const int &,const Symbolic &);
  Symbolic operator+(const Symbolic &,const int &);
  Symbolic operator+(const double &,const Symbolic &);
  Symbolic operator+(const Symbolic &,const double &);
  Symbolic operator++(Symbolic &);
  Symbolic operator++(Symbolic &,int);
  Symbolic operator-(const Symbolic &);
  Symbolic operator-(const Symbolic &,const Symbolic &);
  Symbolic operator-(const int &,const Symbolic &);
  Symbolic operator-(const Symbolic &,const int &);
  Symbolic operator-(const double &,const Symbolic &);
  Symbolic operator-(const Symbolic &,const double &);
  Symbolic operator--(Symbolic &);
  Symbolic operator--(Symbolic &,int);
  Symbolic operator*(const Symbolic &,const Symbolic &);
  Symbolic operator*(const int &,const Symbolic &);
  Symbolic operator*(const Symbolic &,const int &);
  Symbolic operator*(const double &,const Symbolic &);
  Symbolic operator*(const Symbolic &,const double &);
  Symbolic operator/(const Symbolic &,const Symbolic &);
  Symbolic operator/(const int &,const Symbolic &);
  Symbolic operator/(const Symbolic &,const int &);
  Symbolic operator/(const double &,const Symbolic &);
  Symbolic operator/(const Symbolic &,const double &);
  Symbolic operator+=(Symbolic &,const Symbolic &);
  Symbolic operator+=(Symbolic &,const int &);
  Symbolic operator+=(Symbolic &,const double &);
  Symbolic operator-=(Symbolic &,const Symbolic &);
  Symbolic operator-=(Symbolic &,const int &);
  Symbolic operator-=(Symbolic &,const double &);
  Symbolic operator*=(Symbolic &,const Symbolic &);
  Symbolic operator*=(Symbolic &,const int &);
  Symbolic operator*=(Symbolic &,const double &);
  Symbolic operator/=(Symbolic &,const Symbolic &);
  Symbolic operator/=(Symbolic &,const int &);
  Symbolic operator/=(Symbolic &,const double &);
  Equation operator==(const Symbolic &,const Symbolic &);
  Equation operator==(const Symbolic &,int);
  Equation operator==(int,const Symbolic &);
  Equation operator==(const Symbolic &,double);
  Equation operator==(double,const Symbolic &);
  int operator!=(const Symbolic &,const Symbolic &);
  int operator!=(const Symbolic &,int);
  int operator!=(int,const Symbolic &);
  int operator!=(const Symbolic &,double);
  int operator!=(double,const Symbolic &);
  Symbolic sin(const Symbolic &);
  Symbolic cos(const Symbolic &);
  Symbolic tan(const Symbolic &);
  Symbolic cot(const Symbolic &);
  Symbolic sec(const Symbolic &);
  Symbolic csc(const Symbolic &);
  Symbolic sinh(const Symbolic &);
  Symbolic cosh(const Symbolic &);
  Symbolic ln(const Symbolic &);
  Symbolic log(const Symbolic &,const Symbolic &);
  Symbolic pow(const Symbolic &,const Symbolic &);
  Symbolic operator^(const Symbolic &,const Symbolic &);
  Symbolic operator^(const Symbolic &,int);
  Symbolic operator^(int,const Symbolic &);
  Symbolic operator^(const Symbolic &,double);
  Symbolic operator^(double,const Symbolic &);
  Symbolic exp(const Symbolic &);
  Symbolic sqrt(const Symbolic &);
  Symbolic factorial(const Symbolic &);
  Symbolic gamma(const Symbolic &);
  Symbolic df(const Symbolic &,const Symbolic &);
  Symbolic df(const Symbolic &,const Symbolic &,unsigned int);
  Symbolic &rhs(Equations &,const Symbolic &);
  Symbolic &lhs(Equations &,const Symbolic &);
  template<> Symbolic zero(Symbolic);
  template<> Symbolic one(Symbolic);
  Equations operator,(const Equation &,const Equation &);
  Equations operator,(const Equations &,const Equation &);
  Equations operator,(const Equation &,const Equations &);
  std::list<Symbolic> operator,(const Symbolic &,const Symbolic &);
  std::list<Symbolic> operator,(const int &,const Symbolic &);
  std::list<Symbolic> operator,(const double &,const Symbolic &);
  std::list<Symbolic> operator,(const Symbolic &,const int &);
  std::list<Symbolic> operator,(const Symbolic &,const double &);
  std::list<Symbolic> operator,(const std::list<Symbolic> &,const Symbolic &);
  std::list<Symbolic> operator,(const std::list<Symbolic> &,const int &);
  std::list<Symbolic> operator,(const std::list<Symbolic> &,const double &);
  std::list<Symbolic> operator,(const Symbolic &,const std::list<Symbolic> &);
  std::list<Symbolic> operator,(const int &,const std::list<Symbolic> &);
  std::list<Symbolic> operator,(const double &,const std::list<Symbolic> &);
  std::list<std::list<Symbolic> >
    operator,(const std::list<Symbolic> &,const std::list<Symbolic> &);
  std::list<std::list<Symbolic> >
    operator,(const std::list<std::list<Symbolic> > &,const std::list<Symbolic> &);
  std::list<std::list<Symbolic> >
    operator,(const std::list<Symbolic> &,const std::list<std::list<Symbolic> > &);
  Equation operator,(const Symbolic &, const Equation &);
  Equation operator,(const std::list<Symbolic> &, const Equation &);
  ostream &operator<<(ostream &,const Equations &);
  ostream &operator<<(ostream &,const std::list<Symbolic> &);
  Symbolic tr(const Symbolic &);
  Symbolic trace(const Symbolic &);
  Symbolic det(const Symbolic &);
  Symbolic determinant(const Symbolic &);
  Symbolic kron(const Symbolic &,const Symbolic &);
  Symbolic dsum(const Symbolic &,const Symbolic &);
  Symbolic hadamard(const Symbolic &,const Symbolic &);
  void pattern_match_TRUE(PatternMatches &);
  void pattern_match_FALSE(PatternMatches &);
  int pattern_match_AND(Equations&, const Equation&);
  int pattern_match_AND(Equations&, const Equations&);
  void pattern_match_AND(PatternMatches &, const Equation&);
  void pattern_match_AND(PatternMatches &, const Equations&);
  void pattern_match_AND(PatternMatches &, const PatternMatches &);
  void pattern_match_OR(PatternMatches &, const Equation&);
  void pattern_match_OR(PatternMatches &, const Equations&);
  void pattern_match_OR(PatternMatches &, const PatternMatches &);
} //namespace symbolic //DX20200625

//DX20200825 - moved function definitions to symbolicc++.cpp

#endif

