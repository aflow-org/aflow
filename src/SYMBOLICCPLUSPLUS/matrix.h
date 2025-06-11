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


// matrix.h
// Matrix class

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <utility>
#include "identity.h"
#include "vector.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic{ //DX20200625
  // definition of class Matrix
  template <class T> class Matrix
  {
    protected:
      // Data Fields
      int rowNum, colNum;
      Vector<Vector<T> > mat;

    public:
      // Constructors
      Matrix();
      Matrix(int,int);
      Matrix(int,int,const T&);
      Matrix(const Vector<T>&);
      Matrix(const Matrix<T>&);
      ~Matrix();

      // Member Functions
      Vector<T>& operator [] (int);
      const Vector<T>& operator [] (int) const;
      Vector<T>  operator () (int) const;

      Matrix<T> identity();
      Matrix<T> transpose() const;
      Matrix<T> inverse() const;
      T trace() const;
      T determinant() const;

      int rows() const; 
      int cols() const; 
      void resize(int,int);
      void resize(int,int,const T&);
      void fill(const T&);

      // Arithmetic Operators
      const Matrix<T>& operator = (const Matrix<T>&);
      const Matrix<T>& operator = (const T&);

      Matrix<T> operator + () const;
      Matrix<T> operator - () const;
      Matrix<T> operator += (const Matrix<T>&);
      Matrix<T> operator -= (const Matrix<T>&);
      Matrix<T> operator *= (const Matrix<T>&);
      Matrix<T> operator +  (const Matrix<T>&) const;
      Matrix<T> operator -  (const Matrix<T>&) const;
      Matrix<T> operator *  (const Matrix<T>&) const;
      Vector<T> operator *  (const Vector<T>&) const;

      Matrix<T> operator += (const T&);
      Matrix<T> operator -= (const T&);
      Matrix<T> operator *= (const T&);
      Matrix<T> operator /= (const T&);
      Matrix<T> operator +  (const T&) const;
      Matrix<T> operator -  (const T&) const;
      Matrix<T> operator *  (const T&) const;
      Matrix<T> operator /  (const T&) const;

      Vector<T> vec() const;
      Matrix<T> kron(const Matrix<T>&) const;
      Matrix<T> dsum(const Matrix<T>&) const;
      Matrix<T> hadamard(const Matrix<T>&) const;
      std::pair<Matrix<T>, Matrix<T> > LU() const;

      ostream &output(ostream&) const;
      istream &input(istream&);
  };

  //DX20200824 - turned the two lines below into declarations and moved definitions to matrix.cpp
  template <class T> T tr(const Matrix<T> &m);
  template <class T> T det(const Matrix<T> &m);

  template <class T> 
    Matrix<T> hadamard(const Matrix<T> &s,const Matrix<T> &m);

  template <class T> 
    int operator == (const Matrix<T> &m1,const Matrix<T> &m2);

  template <class T> 
    int operator != (const Matrix<T> &m1,const Matrix<T> &m2);


  template <class T> ostream & operator << (ostream &s,const Matrix<T> &m);

  template <class T> istream & operator >> (istream &s,Matrix<T> &m);
} //namespace symbolic //DX20200625
#endif
