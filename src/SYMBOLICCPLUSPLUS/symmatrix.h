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


// symmatrix.h

#ifndef SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX

#include "matrix.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

//DX20200831 [OBSOLETE] #ifdef  SYMBOLIC_FORWARD
//DX20200831 [OBSOLETE] #ifndef SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_FORWARD
//DX20200831 [OBSOLETE] #define SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_FORWARD

namespace symbolic{ //DX20200625
  class SymbolicMatrix;
} //namespace symbolic //DX20200625

//DX20200831 [OBSOLETE] #endif
//DX20200831 [OBSOLETE] #endif

#ifdef  SYMBOLIC_DECLARE
#ifndef SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_DECLARE
#define SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_DECLARE

namespace symbolic{ //DX20200625
  class SymbolicMatrix
    : public CloningSymbolicInterface, public Matrix<Symbolic>
  {
    public: SymbolicMatrix(const SymbolicMatrix&);
            SymbolicMatrix(const Matrix<Symbolic>&);
            SymbolicMatrix(const std::list<std::list<Symbolic> >&);
            SymbolicMatrix(const string&,int,int);
            SymbolicMatrix(const char*,int,int);
            SymbolicMatrix(const Symbolic&,int,int);
            SymbolicMatrix(int,int);
            ~SymbolicMatrix();

            void copy(const SymbolicMatrix&); //DX20210420 - to fix warnings for gcc>10, need explicit declaration
            SymbolicMatrix &operator=(const SymbolicMatrix&); //DX20210420 - to fix warnings for gcc>10, need explicit declaration

            void print(ostream&) const;
            Symbolic subst(const Symbolic&,const Symbolic&,int &n) const;
            Simplified simplify() const;
            int compare(const Symbolic&) const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;
            Symbolic coeff(const Symbolic&) const;
            Expanded expand() const;
            int commute(const Symbolic&) const;
            PatternMatches match(const Symbolic&, const std::list<Symbolic>&) const;
            PatternMatches match_parts(const Symbolic&,
                const std::list<Symbolic>&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };
} //namespace symbolic //DX20200625

#endif
#endif

//DX20200824 - moved function definitions to symmatrix.cpp

#endif
