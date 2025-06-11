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


// matnorm.h
// Norms of Matrices

#ifndef MATNORM_H
#define MATNORM_H

#include <iostream>
#include <math.h>
#include "vector.h"
#include "vecnorm.h"
#include "matrix.h"
#include "identity.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

//DX20200824 - moved function definitions to matnorm.cpp
namespace symbolic{ //DX20200625
  template <class T> T norm1(const Matrix<T> &m);
  template <class T> T normI(const Matrix<T> &m);
  template <class T> T normH(const Matrix<T> &m);
} //namespace symbolic //DX20200625
#endif
