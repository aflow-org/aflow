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


// functions.h

#ifndef SYMBOLIC_CPLUSPLUS_FUNCTIONS

#include <cmath>
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

//DX20200831 [OBSOLETE] #ifdef  SYMBOLIC_FORWARD
//DX20200831 [OBSOLETE] #ifndef SYMBOLIC_CPLUSPLUS_FUNCTIONS_FORWARD
//DX20200831 [OBSOLETE] #define SYMBOLIC_CPLUSPLUS_FUNCTIONS_FORWARD

namespace symbolic{ //DX20200625
  class Sin;
  class Cos;
  class Sinh;
  class Cosh;
  class Log;
  class Power;
  class Derivative;
  class Integral;
  class Rows;
  class Columns;
  class Row;
  class Column;
  class Transpose;
  class Trace;
  class Determinant;
  class Vec;
  class Kronecker;
  class DirectSum;
  class Hadamard;
  class Gamma;
} //namespace symbolic //DX20200625

//DX20200831 [OBSOLETE] #endif
//DX20200831 [OBSOLETE] #endif

#ifdef  SYMBOLIC_DECLARE
#ifndef SYMBOLIC_CPLUSPLUS_FUNCTIONS_DECLARE
#define SYMBOLIC_CPLUSPLUS_FUNCTIONS_DECLARE

namespace symbolic{ //DX20200625
  class Sin: public Symbol
  {
    public: Sin(const Sin&);
            Sin(const Symbolic&);

            Simplified simplify() const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200825
  };

  class Cos: public Symbol
  {
    public: Cos(const Cos&);
            Cos(const Symbolic&);

            Simplified simplify() const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Sinh: public Symbol
  {
    public: Sinh(const Sinh&);
            Sinh(const Symbolic&);

            Simplified simplify() const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Cosh: public Symbol
  {
    public: Cosh(const Cosh&);
            Cosh(const Symbolic&);

            Simplified simplify() const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Log: public Symbol
  {
    public: Log(const Log&);
            Log(const Symbolic&, const Symbolic&);

            void print(ostream&) const;
            Simplified simplify() const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Power: public Symbol
  {
    public: Power(const Power&);
            Power(const Symbolic&,const Symbolic&);

            void print(ostream&) const;
            Simplified simplify() const;
            Expanded expand() const;
            Symbolic subst(const Symbolic &x,const Symbolic &y,int &n) const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;
            PatternMatches match_parts(const Symbolic &s, const std::list<Symbolic>&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Derivative: public Symbol
  {
    public: Derivative(const Derivative&);
            Derivative(const Symbolic&,const Symbolic&);

            Symbolic subst(const Symbolic &x,const Symbolic &y,int &n) const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;
            int compare(const Symbolic&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Integral: public Symbol
  {
    public: Integral(const Integral&);
            Integral(const Symbolic&,const Symbolic&);

            Symbolic subst(const Symbolic &x,const Symbolic &y,int &n) const;
            Symbolic df(const Symbolic&) const;
            Symbolic integrate(const Symbolic&) const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Rows: public Symbol
  {
    public: Rows(const Rows&);
            Rows(const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Columns: public Symbol
  {
    public: Columns(const Columns&);
            Columns(const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Row: public Symbol
  {
    private: int row;
    public:  Row(const Row&);
             Row(const Symbolic&,int);

             Simplified simplify() const;

             //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
             Cloning *clone() const; //DX20200528
  };

  class Column: public Symbol
  {
    private: int column;
    public:  Column(const Column&);
             Column(const Symbolic&,int);

             Simplified simplify() const;

             //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
             Cloning *clone() const; //DX20200528
  };

  class Transpose: public Symbol
  {
    public: Transpose(const Transpose&);
            Transpose(const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Trace: public Symbol
  {
    public: Trace(const Trace&);
            Trace(const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Determinant: public Symbol
  {
    public: Determinant(const Determinant&);
            Determinant(const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Vec: public Symbol
  {
    public: Vec(const Vec&);
            Vec(const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Kronecker: public Symbol
  {
    public: Kronecker(const Kronecker&);
            Kronecker(const Symbolic&,const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class DirectSum: public Symbol
  {
    public: DirectSum(const DirectSum&);
            DirectSum(const Symbolic&,const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Hadamard: public Symbol
  {
    public: Hadamard(const Hadamard&);
            Hadamard(const Symbolic&,const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };

  class Gamma: public Symbol
  {
    public: Gamma(const Gamma&);
            Gamma(const Symbolic&);

            Simplified simplify() const;

            //DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
            Cloning *clone() const; //DX20200528
  };
} //namespace symbolic //DX20200625

#endif
#endif

//DX20200824 - moved function definitions to functions.cpp

#endif
