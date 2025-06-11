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


// multinomial.h

#ifndef _MULTINOMIAL
#define _MULTINOMIAL

#include <cassert>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <utility>  // for pair
#include <vector>
#include "identity.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

// Assumption : T has typecasts defined and can act as a numeric data type

// Multinomial class

namespace symbolic{ //DX20200625
  template <class T>
    class Multinomial
    {
      public:
        Multinomial(): type(number),n(zero(T())) {}
        Multinomial(T);
        Multinomial(string x);
        Multinomial(const Multinomial<T> &p):
          variable(p.variable),type(p.type),n(p.n),u(p.u),m(p.m) {}

        Multinomial<T>& operator=(const T&);
        Multinomial<T>& operator=(const Multinomial<T>&);

        Multinomial<T> operator+() const;
        Multinomial<T> operator-() const;

        Multinomial<T> operator+(const Multinomial<T>&) const;
        Multinomial<T> operator-(const Multinomial<T>&) const;
        Multinomial<T> operator*(const Multinomial<T>&) const;
        Multinomial<T> operator+(const T&) const;
        Multinomial<T> operator-(const T&) const;
        Multinomial<T> operator*(const T&) const;

        Multinomial<T> operator^(unsigned int) const;

        Multinomial<T>& operator+=(const Multinomial<T>&);
        Multinomial<T>& operator-=(const Multinomial<T>&);
        Multinomial<T>& operator*=(const Multinomial<T>&);
        Multinomial<T>& operator+=(const T&);
        Multinomial<T>& operator-=(const T&);
        Multinomial<T>& operator*=(const T&);

        int operator==(const Multinomial<T>&) const;
        int operator!=(const Multinomial<T>&) const;
        int operator==(const T&) const;
        int operator!=(const T&) const;

        Multinomial<T> Diff(const string &) const;
        ostream &output(ostream &) const;

      protected:
        void remove_zeros(void);
        std::pair<Multinomial<T>,Multinomial<T> >
          reconcile(const Multinomial<T>&,const Multinomial<T>&) const;
        vector<string> toarray(void) const;
        string variable;
        enum { number, univariate, multivariate } type;
        T n;                               //number
        std::list<std::pair<T,int> > u;              //univariate
        std::list<std::pair<Multinomial<T>,int> > m; //multivariate
    };


  // additional functions that are not members of the class
  template <class T>
    Multinomial<T> operator+(const T&,const Multinomial<T>&);
  template <class T>
    Multinomial<T> operator-(const T&,const Multinomial<T>&);
  template <class T>
    Multinomial<T> operator*(const T&,const Multinomial<T>&);
  template <class T>
    int operator==(T,const Multinomial<T>&);
  template <class T>
    int operator!=(T,const Multinomial<T>&);

} //namespace symbolic //DX20200625

//DX20200825 - moved function definitions to multinomial.cpp

#endif
