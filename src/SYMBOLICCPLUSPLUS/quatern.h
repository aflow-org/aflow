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


// quatern.h
// Template Class for Quaternions


#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>
#include <cmath>       // for sqrt()
#include "identity.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic{ //DX20200625
  template <class T> class Quaternion
  {
    private:
      // Data Fields      
      T r, i, j, k;

    public:
      // Constructors
      Quaternion();
      Quaternion(T,T,T,T);
      Quaternion(const Quaternion<T>&);
      ~Quaternion();

      // Operators
      const Quaternion<T> &operator = (const Quaternion<T>&);
      Quaternion<T> operator + (const Quaternion<T>&);
      Quaternion<T> operator - (const Quaternion<T>&);
      Quaternion<T> operator - () const;
      Quaternion<T> operator * (const Quaternion<T>&);
      Quaternion<T> operator * (T);
      Quaternion<T> operator / (const Quaternion<T>&);
      Quaternion<T> operator ~ () const;

      // Member Functions
      Quaternion<T> sqr();
      Quaternion<T> conjugate() const;
      Quaternion<T> inverse() const;
      double magnitude() const;

      // Streams
      ostream &print(ostream &) const;
      istream &input(istream &s);
  };

} //namespace symbolic //DX20200625

//DX20200825 - moved function definitions to quatern.cpp

#endif
