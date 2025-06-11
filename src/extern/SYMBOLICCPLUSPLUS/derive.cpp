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

#ifndef SYMBOLIC_CPLUSPLUS_DERIVE_DEFINE //DX20200824 - added for proper compilation
#define SYMBOLIC_CPLUSPLUS_DERIVE_DEFINE //DX20200824 - added for proper compilation

#include "derive.h"

namespace symbolic{ //DX20200824 - symbolic
  template <class T> Derive<T>::Derive() : u(zero(T())),du(one(T())) {}

  template <class T> Derive<T>::Derive(const T &v) : u(v),du(zero(T())) {}

  template <class T> Derive<T>::Derive(const T &v,const T &dv) : u(v),du(dv) {}

  template <class T> Derive<T>::Derive(const Derive<T> &r) : u(r.u),du(r.du) {}

  template <class T> void Derive<T>::set(const T v) { u = v; }

  template <class T> Derive<T> Derive<T>::operator - () const
  { return Derive<T>(-u,-du); }

  template <class T> Derive<T> Derive<T>::operator += (const Derive<T> &r)
  { return *this = *this + r; }

  template <class T> Derive<T> Derive<T>::operator -= (const Derive<T> &r)
  { return *this = *this - r; }

  template <class T> Derive<T> Derive<T>::operator *= (const Derive<T> &r)
  { return *this = *this * r; }

  template <class T> Derive<T> Derive<T>::operator /= (const Derive<T> &r)
  { return *this = *this / r; }

  template <class T> Derive<T> Derive<T>::operator += (const T &r)
  { return *this = *this + r; }

  template <class T> Derive<T> Derive<T>::operator -= (const T &r)
  { return *this = *this - r; }

  template <class T> Derive<T> Derive<T>::operator *= (const T &r)
  { return *this = *this * r; }

  template <class T> Derive<T> Derive<T>::operator /= (const T &r)
  { return *this = *this / r; }

  template <class T>
    Derive<T> Derive<T>::operator + (const Derive<T> &y) const
    { return Derive<T>(u+y.u,du+y.du); }

  template <class T>
    Derive<T> Derive<T>::operator - (const Derive<T> &y) const
    { return Derive<T>(u-y.u,du-y.du); }

  template <class T>
    Derive<T> Derive<T>::operator * (const Derive<T> &y) const
    { return Derive<T>(u*y.u,y.u*du+u*y.du); }

  template <class T>
    Derive<T> Derive<T>::operator / (const Derive<T> &y) const
    { return Derive<T>(u/y.u,(y.u*du-u*y.du)/(y.u*y.u)); }

  template <class T>
    Derive<T> Derive<T>::operator + (const T &y) const
    { return *this + Derive<T>(y); }

  template <class T>
    Derive<T> Derive<T>::operator - (const T &y) const
    { return *this - Derive<T>(y); }

  template <class T>
    Derive<T> Derive<T>::operator * (const T &y) const
    { return *this * Derive<T>(y); }

  template <class T>
    Derive<T> Derive<T>::operator / (const T &y) const
    { return *this / Derive<T>(y); }

  template <class T>
    Derive<T> operator + (const T &y,const Derive<T> &d)
    { return Derive<T>(y)+d; }

  template <class T>
    Derive<T> operator - (const T &y,const Derive<T> &d)
    { return Derive<T>(y)-d; }

  template <class T>
    Derive<T> operator * (const T &y,const Derive<T> &d)
    { return Derive<T>(y)*d; }

  template <class T>
    Derive<T> operator / (const T &y,const Derive<T> &d)
    { return Derive<T>(y)/d; }

  template <class T>
    Derive<T> Derive<T>::exp() const { return Derive<T>(exp(u),du*exp(u)); }

  template <class T>
    Derive<T> exp(const Derive<T> &x) { return x.exp(); }

  template <class T>
    Derive<T> Derive<T>::sin() const { return Derive<T>(sin(u),du*cos(u)); }

  template <class T>
    Derive<T> sin(const Derive<T> &x) { return x.sin(); }

  template <class T>
    Derive<T> Derive<T>::cos() const { return Derive<T>(cos(u),-du*sin(u)); }

  template <class T>
    Derive<T> cos(const Derive<T> &x) { return x.cos(); }

  template <class T> T Derive<T>::df() const { return du; }

  template <class T> T df(const Derive<T> &x) { return x.df(); }

  template <class T>
    ostream &Derive<T>::output(ostream &s) const { return s << u; }

  template <class T>
    ostream &operator << (ostream &s,const Derive<T> &r)
    { return r.output(s); }
} //DX20200824 - symbolic namespace

#endif //DX20200824 - added for proper compilation
