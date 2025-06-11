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

#ifndef SYMBOLIC_CPLUSPLUS_VECNORM_DEFINE //DX20200824 - added for proper compilation
#define SYMBOLIC_CPLUSPLUS_VECNORM_DEFINE //DX20200824 - added for proper compilation
#include "vecnorm.h"

namespace symbolic{ //DX20200625
  template <class T> T norm1(const Vector<T> &v)
  {
    T result(0);
    for(int i=0;i<int(v.size());i++) result = result + abs(v[i]);
    return result;
  }

} //namespace symbolic //DX20200625

namespace symbolic{ //DX20200625
  double norm1(const Vector<double> &v)
  {
    double result(0);
    for(int i=0;i<int(v.size());i++) result = result + fabs(v[i]);
    return result;
  }
} //namespace symbolic //DX20200625

namespace symbolic{ //DX20200625
  template <class T> double norm2(const Vector<T> &v)
  {
    T result(0);
    for(int i=0;i<int(v.size());i++) result = result + v[i]*v[i];
    return std::sqrt(double(result)); //DX20200728 - need to specify std otherwise it is ambiguous
  }

  template <class T> T normI(const Vector<T> &v)
  {
    T maxItem(abs(v[0])), temp;
    for(int i=1;i<int(v.size());i++)
    {
      temp = abs(v[i]);
      if(temp > maxItem) maxItem = temp;
    }
    return maxItem;
  }

} //namespace symbolic //DX20200625

namespace symbolic{ //DX20200625
  double normI(const Vector<double> &v)
  {
    double maxItem(fabs(v[0])), temp;
    for(int i=1;i<int(v.size());i++)
    {
      temp = fabs(v[i]);
      if(temp > maxItem) maxItem = temp;
    }
    return maxItem;
  }
} //namespace symbolic //DX20200625

namespace symbolic{ //DX20200625
  template <class T> Vector<T> normalize(const Vector<T> &v)
  {
    Vector<T> result(v.size());
    double length = norm2(v);
    for(int i=0;i<int(v.size());i++) result[i] = v[i]/length;
    return result;
  }
} //namespace symbolic //DX20200625
#endif //DX20200824 - added for proper compilation
