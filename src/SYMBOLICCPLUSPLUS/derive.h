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


// derive.h
// The Derivation Class

#ifndef DERIVE_H
#define DERIVE_H

#include <iostream>
#include <cmath>
#include "identity.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic{ //DX20200625
  template <class T> class Derive
  {
    private:
      // Data Field
      T u, du;

      // Private Constructor
      Derive(const T&,const T&);

    public:
      // Constructors
      Derive();
      Derive(const T&);
      Derive(const Derive<T>&);

      // Member Function
      void set(const T);

      // Arithmetic Operators
      Derive<T> operator - () const;
      Derive<T> operator += (const Derive<T>&);
      Derive<T> operator -= (const Derive<T>&);
      Derive<T> operator *= (const Derive<T>&);
      Derive<T> operator /= (const Derive<T>&);
      Derive<T> operator += (const T&);
      Derive<T> operator -= (const T&);
      Derive<T> operator *= (const T&);
      Derive<T> operator /= (const T&);

      Derive<T> operator + (const Derive<T>&) const;
      Derive<T> operator - (const Derive<T>&) const;
      Derive<T> operator * (const Derive<T>&) const;
      Derive<T> operator / (const Derive<T>&) const;
      Derive<T> operator + (const T&) const;
      Derive<T> operator - (const T&) const;
      Derive<T> operator * (const T&) const;
      Derive<T> operator / (const T&) const;
      Derive<T> exp() const;
      Derive<T> sin() const;
      Derive<T> cos() const;
      T df() const;
      ostream &output(ostream&) const;
  };

} //namespace symbolic //DX20200625

#endif
