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

// symbolic.h

#ifndef SYMBOLIC_CPLUSPLUS_SYMBOLIC

#include <iostream>
#include <list>
#include <typeinfo>
// using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic { // DX20200625
  class CloningSymbolicInterface;
  class Expanded;
  class Simplified;
  class Symbolic;
  class SymbolicInterface;
  class SymbolicProxy;
} // namespace symbolic

#ifdef SYMBOLIC_DECLARE
#ifndef SYMBOLIC_CPLUSPLUS_SYMBOLIC_DECLARE
#define SYMBOLIC_CPLUSPLUS_SYMBOLIC_DECLARE

namespace symbolic { // DX20200625
  class SymbolicInterface {
  public:
    int simplified, expanded;
    SymbolicInterface();
    SymbolicInterface(const SymbolicInterface &);
    virtual ~SymbolicInterface();

    virtual void print(std::ostream &) const = 0;
    [[nodiscard]] virtual const std::type_info &type() const;
    virtual Symbolic subst(const Symbolic &, const Symbolic &, int &n) const = 0;
    [[nodiscard]] virtual Simplified simplify() const = 0;
    [[nodiscard]] virtual int compare(const Symbolic &) const = 0;
    [[nodiscard]] virtual Symbolic df(const Symbolic &) const = 0;
    [[nodiscard]] virtual Symbolic integrate(const Symbolic &) const = 0;
    [[nodiscard]] virtual Symbolic coeff(const Symbolic &) const = 0;
    [[nodiscard]] virtual Expanded expand() const = 0;
    [[nodiscard]] virtual int commute(const Symbolic &) const = 0;
            // match *this (as a pattern) against an expression
    [[nodiscard]] virtual PatternMatches match(const Symbolic &, const std::list<Symbolic> &) const = 0;
            // match parts of *this (as an expression) against a pattern
    [[nodiscard]] virtual PatternMatches match_parts(const Symbolic &, const std::list<Symbolic> &) const = 0;
  };

  class CloningSymbolicInterface : public SymbolicInterface, public Cloning {
  public:
    CloningSymbolicInterface();
    CloningSymbolicInterface(const CloningSymbolicInterface &);
    void copy(const CloningSymbolicInterface &c); // DX20211120 - to fix warnings for gcc>10 and gcc=4, need explicit declaration
    CloningSymbolicInterface &operator=(const CloningSymbolicInterface &); // DX20210420 - to fix warnings for gcc>10, need explicit declaration
  };

  class SymbolicProxy : public SymbolicInterface, public CastPtr<CloningSymbolicInterface> {
  public:
    SymbolicProxy(const CloningSymbolicInterface &);
    SymbolicProxy(const SymbolicProxy &);
    SymbolicProxy(const Number<void> &);
    SymbolicProxy();

    void print(std::ostream &) const;
    [[nodiscard]] const std::type_info &type() const;
    Symbolic subst(const Symbolic &, const Symbolic &, int &n) const;
    [[nodiscard]] Simplified simplify() const;
    [[nodiscard]] int compare(const Symbolic &) const;
    [[nodiscard]] Symbolic df(const Symbolic &) const;
    [[nodiscard]] Symbolic integrate(const Symbolic &) const;
    [[nodiscard]] Symbolic coeff(const Symbolic &) const;
    [[nodiscard]] Expanded expand() const;
    [[nodiscard]] int commute(const Symbolic &) const;
    [[nodiscard]] PatternMatches match(const Symbolic &, const std::list<Symbolic> &) const;
    [[nodiscard]] PatternMatches match_parts(const Symbolic &, const std::list<Symbolic> &) const;

    SymbolicProxy &operator=(const CloningSymbolicInterface &);
    SymbolicProxy &operator=(const SymbolicProxy &);
  };

  class Simplified : public SymbolicProxy {
  public:
    Simplified(const CloningSymbolicInterface &);
    Simplified(const SymbolicProxy &);
    Simplified(const Number<void> &);
  };

  class Expanded : public SymbolicProxy {
  public:
    Expanded(const CloningSymbolicInterface &);
    Expanded(const SymbolicProxy &);
    Expanded(const Number<void> &);
  };

  class Symbolic : public SymbolicProxy {
  public:
    static int auto_expand;
    static int subst_count;

    Symbolic();
    Symbolic(const Symbolic &);
    Symbolic(const CloningSymbolicInterface &);
    Symbolic(const SymbolicProxy &);
    Symbolic(const Number<void> &);
    Symbolic(const int &);
    Symbolic(const double &);
    Symbolic(const std::string &);
    Symbolic(const char *);
    Symbolic(const std::string &, int);
    Symbolic(const char *, int);
    Symbolic(const Symbolic &, int);
    Symbolic(const std::string &, int, int);
    Symbolic(const char *, int, int);
    Symbolic(const Symbolic &, int, int);
    Symbolic(const std::list<Symbolic> &);
    Symbolic(const std::list<std::list<Symbolic>> &);
    ~Symbolic();

    SymbolicProxy &operator=(const CloningSymbolicInterface &);
    SymbolicProxy &operator=(const SymbolicProxy &);
    SymbolicProxy &operator=(const Symbolic &); // DX20200908 - missing this assignment operator in main SYMBOLICC++ source
    SymbolicProxy &operator=(const int &);
    SymbolicProxy &operator=(const double &);
    SymbolicProxy &operator=(const std::string &);
    SymbolicProxy &operator=(const char *);
    SymbolicProxy &operator=(const std::list<Symbolic> &);
    SymbolicProxy &operator=(const std::list<std::list<Symbolic>> &);
    Symbolic operator[](const Equation &) const;
    Symbolic operator[](const Equations &) const;
    Symbolic operator[](const Symbolic &) const;
    Symbolic operator[](const std::list<Symbolic> &) const;
    Symbolic &operator()(int);
    Symbolic &operator()(int, int);
    const Symbolic &operator()(int) const;
    const Symbolic &operator()(int, int) const;
    Symbolic subst(const Symbolic &, const Symbolic &, int &n = subst_count) const;
    Symbolic subst(const Symbolic &, const int &, int &n = subst_count) const;
    Symbolic subst(const Symbolic &, const double &, int &n = subst_count) const;
    Symbolic subst(const Equation &, int &n = subst_count) const;
    Symbolic subst(const Equations &, int &n = subst_count) const;
    Symbolic subst_all(const Symbolic &, const Symbolic &, int &n = subst_count) const;
    Symbolic subst_all(const Equation &, int &n = subst_count) const;
    Symbolic subst_all(const Equations &, int &n = subst_count) const;
    [[nodiscard]] Symbolic coeff(const Symbolic &) const;
    [[nodiscard]] Symbolic coeff(const Symbolic &, int) const;
    [[nodiscard]] Symbolic coeff(const int &) const;
    [[nodiscard]] Symbolic coeff(const double &) const;

    [[nodiscard]] Symbolic commutative(int) const;
    Symbolic operator~() const;
    operator int() const;
    operator double() const;

    Symbolic operator|(const Symbolic &) const;
    Symbolic operator%(const Symbolic &) const;

    [[nodiscard]] int rows() const;
    [[nodiscard]] int columns() const;
    [[nodiscard]] Symbolic row(int) const;
    [[nodiscard]] Symbolic column(int) const;
    [[nodiscard]] Symbolic identity() const;
    [[nodiscard]] Symbolic transpose() const;
    [[nodiscard]] Symbolic trace() const;
    [[nodiscard]] Symbolic determinant() const;
    [[nodiscard]] Symbolic vec() const;
    [[nodiscard]] Symbolic kron(const Symbolic &) const;
    [[nodiscard]] Symbolic dsum(const Symbolic &) const;
    [[nodiscard]] Symbolic hadamard(const Symbolic &) const;
    [[nodiscard]] Symbolic inverse() const;
  };
} // namespace symbolic

#endif
#endif

#endif
