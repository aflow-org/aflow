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

//DX20200831 - created this file and moved function definitions here and removed unnecessary macros

#ifndef SYMBOLIC_CPLUSPLUS_IDENTITY_DEFINE //DX20200824 - added for proper compilation
#define SYMBOLIC_CPLUSPLUS_IDENTITY_DEFINE //DX20200824 - added for proper compilation

#include "identity.h"

namespace symbolic{ //DX20200625
  template <class T> T zero(T) { static const T z = T() - T(); return z; }

  template <class T> T one(T)
  {
    static T z(zero(T()));
    static const T o = ++z;
    return o;
  }
}

namespace symbolic{ //DX20200625
  template <> char zero(char) { return (char) 0; }
  template <> char one(char) { return (char) 1; }

  template <> short zero(short) { return (short) 0; }
  template <> short one(short) { return (short) 1; }

  template <> int zero(int) { return (int) 0; }
  template <> int one(int) { return (int) 1; }

  template <> long zero(long) { return (long) 0; }
  template <> long one(long) { return (long) 1; }

  template <> float zero(float) { return (float) 0.0; }
  template <> float one(float) { return (float) 1.0; }

  template <> double zero(double) { return (double) 0.0; }
  template <> double one(double) { return (double) 1.0; }
} //namespace symbolic //DX20200625


namespace symbolic{ //DX20200625
  template <class T>
    std::complex<T> zero(std::complex<T>) { return std::complex<T>(zero(T())); }

  template <class T>
    std::complex<T> one(std::complex<T>) { return std::complex<T>(one(T())); }
} //namespace symbolic //DX20200625

#endif //DX20200824 - added for proper compilation
