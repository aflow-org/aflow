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


// vecnorm.h
// Norms of Vectors

#ifndef MVECNORM_H
#define MVECNORM_H

#include <iostream>
#include <math.h>
#include "vector.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

//DX20200825 - moved all function definitions to vecnorm.cpp

namespace symbolic{ //DX20200625
  template <class T> T norm1(const Vector<T> &v);
  double norm1(const Vector<double> &v);
  template <class T> double norm2(const Vector<T> &v);
  template <class T> T normI(const Vector<T> &v);
  double normI(const Vector<double> &v);
  template <class T> Vector<T> normalize(const Vector<T> &v);
} //namespace symbolic //DX20200625
#endif
