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


// array.h
// The Array Class

#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <cassert>
#include <cstdarg>
#include <vector>
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic{ //DX20200625
  //DX20200824 - turned the line below into a declaration and moved the definition into array.cpp
  template <int d> vector<int> dimensions(int d1, ...);

  template <class T,int d>
    class Array : public vector<Array<T,d-1> >
  {
    public:
      // Constructors
      Array();
      Array(int,...);
      Array(vector<int>);
      Array(vector<int>,const T&);
      Array(const Array<T,d>&);
      ~Array();

      // Member Functions
      void resize(int,...);
      void resize(vector<int>);
      void resize(vector<int>,const T&);

      // Arithmetic Operators
      const Array<T,d> & operator = (const T&);
      Array<T,d> operator *= (const T&);
      Array<T,d> operator += (const Array<T,d>&);
      Array<T,d> operator -= (const Array<T,d>&);
      Array<T,d> operator +  (const Array<T,d>&);
      Array<T,d> operator -  (const Array<T,d>&);
  };    // end declaration class Array<T,d>
} //namespace symbolic //DX20200625

//DX20200824 - moved definitions below into array.cpp

#endif
