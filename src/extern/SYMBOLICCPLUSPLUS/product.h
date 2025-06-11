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

// product.h

#ifndef SYMBOLIC_CPLUSPLUS_PRODUCT

#include <algorithm>
#include <list>
#include <vector>
// using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic { // DX20200625
  class Product;
} // namespace symbolic

#ifdef SYMBOLIC_DECLARE
#ifndef SYMBOLIC_CPLUSPLUS_PRODUCT_DECLARE
#define SYMBOLIC_CPLUSPLUS_PRODUCT_DECLARE

namespace symbolic { // DX20200625
  class Product : public CloningSymbolicInterface {
  public:
    std::list<Symbolic> factors;
    Product();
    Product(const Product &);
    Product(const Symbolic &, const Symbolic &);
    ~Product();

    Product &operator=(const Product &);
    [[nodiscard]] int printsNegative() const;

    void print(std::ostream &) const;
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

            // DX20200825 - Cloning *clone() const { return Cloning::clone(*this); }
    [[nodiscard]] Cloning *clone() const; // DX20200528
  };
} // namespace symbolic

#endif
#endif

// DX20200824 - moved function definitions to product.cpp

#endif
