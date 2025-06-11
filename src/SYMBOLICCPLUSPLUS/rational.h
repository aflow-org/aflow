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


// rational.h
// Rational Numbers Class
#ifndef RATIONAL_H
#define RATIONAL_H

#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <limits>
#include <vector>
#include <ctype.h>
#include "identity.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic{ //DX20200625
  template <class T>
    class Rational
    {
      public:
        enum output_type { decimal, fraction };

      private:
        // Data Fields : Numerator and Denominator
        T p,q;

        // Private member function
        T gcd(T,T);

        // output type : decimal or fraction
        static output_type output_preference;

      public:
        // Constructors and Destructor
        Rational();
        Rational(T);
        Rational(T,T);
        Rational(const Rational<T>&);
        Rational(const string&);
        Rational(const double&);
        ~Rational();

        // Conversion operator
        operator double () const;

        // Member functions
        T num() const;            // numerator of r
        T den() const;            // denominator of r
        Rational<T> frac() const; // fractional part of r
        void normalize();         // normalize the rational number

        // Arithmetic operators and Relational operators
        const Rational<T> &operator = (const Rational<T>&);
        Rational<T> operator - () const;
        const Rational<T> &operator += (const Rational<T>&);
        const Rational<T> &operator -= (const Rational<T>&);
        const Rational<T> &operator *= (const Rational<T>&);
        const Rational<T> &operator /= (const Rational<T>&);

        Rational<T> operator + (const Rational<T>&) const;
        Rational<T> operator - (const Rational<T>&) const;
        Rational<T> operator * (const Rational<T>&) const;
        Rational<T> operator / (const Rational<T>&) const;

        Rational<T> operator ++ ();
        Rational<T> operator ++ (int);
        Rational<T> operator -- ();
        Rational<T> operator -- (int);

        int operator == (const Rational<T>&) const;
        int operator != (const Rational<T>&) const;
        int operator >  (const Rational<T>&) const;
        int operator <  (const Rational<T>&) const;
        int operator >= (const Rational<T>&) const;
        int operator <= (const Rational<T>&) const;

        // I/O stream functions
        ostream &output(ostream&) const;
        istream &input(istream&);

        static void output_style(output_type);
    };

  template <class T>
    typename Rational<T>::output_type
    Rational<T>::output_preference = Rational<T>::fraction;

} //namespace symbolic //DX20200625

#include "verylong.h"

namespace symbolic{ //DX20200625
  template <> Rational<Verylong>::operator double() const;
} //namespace symbolic //DX20200625

//DX20200825 - moved function definitions to rational.cpp

#endif
