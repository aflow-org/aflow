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


#ifndef IDENTITY_H
#define IDENTITY_H

#include <iostream>
#include <cstdlib>
#include <complex>
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic{ //DX20200625
  template <class T> T zero(T);
  template <class T> T one(T);

  template <> char zero(char);
  template <> char one(char);
  template <> short zero(short);
  template <> short one(short);
  template <> int zero(int);
  template <> int one(int);
  template <> long zero(long);
  template <> long one(long);
  template <> float zero(float);
  template <> float one(float);
  template <> double zero(double);
  template <> double one(double);
} //namespace symbolic //DX20200625


#endif
