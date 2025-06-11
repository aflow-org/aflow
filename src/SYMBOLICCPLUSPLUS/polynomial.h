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


// polynomial.h

#ifndef _POLYNOMIAL
#define _POLYNOMIAL

#include <cassert>
#include <iostream>
#include <list>
#include <string>
#include <utility>  // for pair
#include "identity.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

//Polynomial class

namespace symbolic{ //DX20200625
  template <class T>
    class Polynomial
    {
      public:
        static int Karatsuba, Newton;

        Polynomial() {}
        Polynomial(const T&);
        Polynomial(string x);
        Polynomial(const Polynomial<T> &p)
          : variable(p.variable), terms(p.terms) {}

        Polynomial<T>& operator=(const T&);
        Polynomial<T>& operator=(const Polynomial<T>&);

        Polynomial<T> operator+() const;
        Polynomial<T> operator-() const;

        Polynomial<T> operator+(const Polynomial<T>&) const;
        Polynomial<T> operator-(const Polynomial<T>&) const;
        Polynomial<T> operator*(const Polynomial<T>&) const;
        Polynomial<T> operator/(const Polynomial<T>&) const;
        Polynomial<T> operator%(const Polynomial<T>&) const;

        Polynomial<T> operator^(unsigned int) const;

        Polynomial<T> operator+(const T&) const;
        Polynomial<T> operator-(const T&) const;
        Polynomial<T> operator*(const T&) const;
        Polynomial<T> operator/(const T&) const;
        Polynomial<T> operator%(const T&) const;

        Polynomial<T>& operator+=(const Polynomial<T>&);
        Polynomial<T>& operator-=(const Polynomial<T>&);
        Polynomial<T>& operator*=(const Polynomial<T>&);
        Polynomial<T>& operator/=(const Polynomial<T>&);
        Polynomial<T>& operator%=(const Polynomial<T>&);

        Polynomial<T>& operator+=(const T&);
        Polynomial<T>& operator-=(const T&);
        Polynomial<T>& operator*=(const T&);
        Polynomial<T>& operator/=(const T&);
        Polynomial<T>& operator%=(const T&);

        Polynomial<T> karatsuba(const Polynomial<T>&,
            const Polynomial<T>&,int) const;
        Polynomial<T> newton(const Polynomial<T>&,const Polynomial<T>&) const;
        Polynomial<T> gcd(const Polynomial<T>&,const Polynomial<T>&) const;
        Polynomial<T> Diff(const string &) const;
        Polynomial<T> Int(const string &) const;
        Polynomial<T> reverse() const;
        std::list<Polynomial<T> > squarefree() const;

        int operator==(const Polynomial<T>&) const;
        int operator!=(const Polynomial<T>&) const;

        int operator==(const T&) const;
        int operator!=(const T&) const;

        T operator()(const T&) const;
        ostream &output1(ostream &) const;
        ostream &output2(ostream &) const;
        int constant() const;
      protected:
        void remove_zeros(void);
        string variable;
        std::list<std::pair<T,int> > terms;
    };

  // multiply using Karatsuba algorithm
  template <class T> int Polynomial<T>::Karatsuba = 1;
  // inversion by Newton iteration for division
  template <class T> int Polynomial<T>::Newton = 1;

  // additional functions that are not members of the class
  template <class T>
    Polynomial<T> operator+(const T&,const Polynomial<T>&);
  template <class T>
    Polynomial<T> operator-(const T&,const Polynomial<T>&);
  template <class T>
    Polynomial<T> operator*(const T&,const Polynomial<T>&);
  template <class T>
    Polynomial<T> operator/(const T&,const Polynomial<T>&);
  template <class T>
    Polynomial<T> operator%(const T&,const Polynomial<T>&);
  template <class T>
    int operator==(T,const Polynomial<T>&);
  template <class T>
    int operator!=(T,const Polynomial<T>&);


} //namespace symbolic //DX20200625

//DX20200825 - moved function definitions to polynomial.cpp

#endif
