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

// number.h

#ifndef SYMBOLIC_CPLUSPLUS_NUMBER

#include <cmath>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <utility>

#include "rational.h"
#include "verylong.h"
// using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic { // DX20200625
  class Numeric;
  template <class T> class Number;
} // namespace symbolic

#ifdef SYMBOLIC_DECLARE
#ifndef SYMBOLIC_CPLUSPLUS_NUMBER_DECLARE
#define SYMBOLIC_CPLUSPLUS_NUMBER_DECLARE

namespace symbolic { // DX20200625
  class Numeric : public CloningSymbolicInterface {
  public:
    Numeric();
    Numeric(const Numeric &);
    [[nodiscard]] virtual const std::type_info &numerictype() const = 0;
    [[nodiscard]] virtual Number<void> add(const Numeric &) const = 0;
    [[nodiscard]] virtual Number<void> mul(const Numeric &) const = 0;
    [[nodiscard]] virtual Number<void> div(const Numeric &) const = 0;
    [[nodiscard]] virtual Number<void> mod(const Numeric &) const = 0;
    [[nodiscard]] virtual int isZero() const = 0;
    [[nodiscard]] virtual int isOne() const = 0;
    [[nodiscard]] virtual int isNegative() const = 0;
    [[nodiscard]] virtual int cmp(const Numeric &) const = 0;
    static std::pair<Number<void>, Number<void>> match(const Numeric &, const Numeric &);
    Symbolic subst(const Symbolic &, const Symbolic &, int &n) const;
    [[nodiscard]] int compare(const Symbolic &) const;
    [[nodiscard]] Symbolic df(const Symbolic &) const;
    [[nodiscard]] Symbolic integrate(const Symbolic &) const;
    [[nodiscard]] Symbolic coeff(const Symbolic &) const;
    [[nodiscard]] Expanded expand() const;
    [[nodiscard]] int commute(const Symbolic &) const;
    [[nodiscard]] PatternMatches match(const Symbolic &, const std::list<Symbolic> &) const;
    [[nodiscard]] PatternMatches match_parts(const Symbolic &, const std::list<Symbolic> &) const;
  };

  template <class T> class Number : public Numeric {
  public:
    T n;
    Number();
    Number(const Number &);
    Number(const T &);
    ~Number();

    Number &operator=(const Number &);
    Number &operator=(const T &);

    void print(std::ostream &) const;
    [[nodiscard]] const std::type_info &type() const;
    [[nodiscard]] const std::type_info &numerictype() const;
    [[nodiscard]] Simplified simplify() const;

    [[nodiscard]] Number<void> add(const Numeric &) const;
    [[nodiscard]] Number<void> mul(const Numeric &) const;
    [[nodiscard]] Number<void> div(const Numeric &) const;
    [[nodiscard]] Number<void> mod(const Numeric &) const;
    [[nodiscard]] int isZero() const;
    [[nodiscard]] int isOne() const;
    [[nodiscard]] int isNegative() const;
    [[nodiscard]] int cmp(const Numeric &) const;

            // DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
    [[nodiscard]] Cloning *clone() const; // DX20200528
  };

  template <> class Number<void> : public CastPtr<Numeric> {
  public:
    Number();
    Number(const Number &);
    Number(const Numeric &);
    Number(const Symbolic &);
    ~Number();

    [[nodiscard]] const std::type_info &numerictype() const;
    [[nodiscard]] int isZero() const;
    [[nodiscard]] int isOne() const;
    [[nodiscard]] int isNegative() const;
    static std::pair<Number<void>, Number<void>> match(const Numeric &, const Numeric &);
    static std::pair<Number<void>, Number<void>> match(const Number<void> &, const Number<void> &);

    void copy(const Numeric &); // DX20210420 - to fix warnings for gcc>10, need explicit declaration
    Number<void> operator+(const Numeric &) const;
    Number<void> operator-(const Numeric &) const;
    Number<void> operator*(const Numeric &) const;
    Number<void> operator/(const Numeric &) const;
    Number<void> operator%(const Numeric &) const;
    Number<void> &operator=(const Numeric &); // DX20210420 - to fix warnings for gcc>10, need explicit declaration
    Number<void> &operator+=(const Numeric &);
    Number<void> &operator*=(const Numeric &);
    Number<void> &operator/=(const Numeric &);
    Number<void> &operator%=(const Numeric &);
    int operator==(const Numeric &) const;
    int operator<(const Numeric &) const;
    int operator>(const Numeric &) const;
    int operator<=(const Numeric &) const;
    int operator>=(const Numeric &) const;

    void copy(const Number<void> &); // DX20210420 - to fix warnings for gcc>10, need explicit declaration
    Number<void> operator+(const Number<void> &) const;
    Number<void> operator-(const Number<void> &) const;
    Number<void> operator*(const Number<void> &) const;
    Number<void> operator/(const Number<void> &) const;
    Number<void> operator%(const Number<void> &) const;
    Number<void> &operator=(const Number<void> &); // DX20210420 - to fix warnings for gcc>10, need explicit declaration
    Number<void> &operator+=(const Number<void> &);
    Number<void> &operator*=(const Number<void> &);
    Number<void> &operator/=(const Number<void> &);
    Number<void> &operator%=(const Number<void> &);
    int operator==(const Number<void> &) const;
    int operator<(const Number<void> &) const;
    int operator>(const Number<void> &) const;
    int operator<=(const Number<void> &) const;
    int operator>=(const Number<void> &) const;

    Numeric &operator*() const;
  };

  Number<void> operator+(const Numeric &, const Number<void> &);
  Number<void> operator-(const Numeric &, const Number<void> &);
  Number<void> operator*(const Numeric &, const Number<void> &);
  Number<void> operator/(const Numeric &, const Number<void> &);
  Number<void> operator%(const Numeric &, const Number<void> &);

  // Template specialization for Rational<Number<void> > and identities
  template <> Rational<Number<void>>::operator double() const;
  template <> Number<void> zero(Number<void>);
  template <> Number<void> one(Number<void>);
} // namespace symbolic
#endif
#endif

// DX20200824 - moved function definitions to number.cpp

#endif
