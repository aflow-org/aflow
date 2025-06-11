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


// vector.h
// Vector class

#ifndef MVECTOR_H
#define MVECTOR_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include "identity.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic{ //DX20200625
  // definition of class Vector
  template <class T> class Vector: public vector<T>
  {
    public:
      // Constructors
      Vector();
      Vector(int);
      Vector(int,const T&);
      Vector(const Vector<T>&);
      ~Vector();

      // Member Functions
      T& operator [] (int);
      const T& operator [] (int) const;
      void reset(int);
      void reset(int,const T&);

      // Arithmetic Operators
      const Vector<T>& operator = (const T&);
      Vector<T> operator + () const;
      Vector<T> operator - () const;

      void copy(const Vector<T>&); //DX20210420 - to fix warnings for gcc>10, need explicit declaration
      Vector<T> &operator = (const Vector<T>&); //DX20210420 - to fix warnings for gcc>10, need explicit declaration
      Vector<T> operator += (const Vector<T>&);
      Vector<T> operator -= (const Vector<T>&);
      Vector<T> operator *= (const Vector<T>&);
      Vector<T> operator /= (const Vector<T>&);
      Vector<T> operator +  (const Vector<T>&) const;
      Vector<T> operator -  (const Vector<T>&) const;
      Vector<T> operator *  (const Vector<T>&) const;
      Vector<T> operator /  (const Vector<T>&) const;

      Vector<T> operator += (const T&);
      Vector<T> operator -= (const T&);
      Vector<T> operator *= (const T&);
      Vector<T> operator /= (const T&);
      Vector<T> operator +  (const T&) const;
      Vector<T> operator -  (const T&) const;
      Vector<T> operator *  (const T&) const;
      Vector<T> operator /  (const T&) const;

      T operator | (const Vector<T>&) const; // Dot product / Inner product
      Vector<T> operator % (const Vector<T>&) const; // Cross product

      ostream &output(ostream&) const;
      istream &input(istream&);
  };

  //DX20200825 - moved function definitions to vector.cpp

} //namespace symbolic //DX20200625
#endif
