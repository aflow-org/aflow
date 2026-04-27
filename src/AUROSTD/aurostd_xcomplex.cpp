// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#include "aurostd_xcomplex.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "aurostd_automatic_template.h"
#include "aurostd_defs.h"
#include "aurostd_xscalar.h"

#ifndef XXEND
#define XXEND 1
#endif

using std::ifstream;
using std::iostream;
using std::istringstream;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::stringstream;

// ----------------------------------------------------------------------------
// --------------------------------------------------------------- constructors
namespace aurostd {

  // constructors
  //   template<class utype>
  //   xcomplex<utype>::xcomplex(void) {
  //     re=utype(0);
  //     im=utype(0);
  //   }

  // destructor
  template <class utype> xcomplex<utype>::~xcomplex() {
    free();
  }
#define AST_TEMPLATE(utype) template xcomplex<utype>::~xcomplex();
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  template <class utype> void xcomplex<utype>::free() {
      // nothing to be done
  }
#define AST_TEMPLATE(utype) template void xcomplex<utype>::free();
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // copy
  template <class utype> xcomplex<utype>::xcomplex(const xcomplex<utype>& b) {
    *this = b;
  }
#define AST_TEMPLATE(utype) template xcomplex<utype>::xcomplex(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  template <class utype> const xcomplex<utype>& xcomplex<utype>::operator=(const xcomplex<utype>& r) {        // copy
    this->re = r.re;
    this->im = r.im;
    return *this;
  }
#define AST_TEMPLATE(utype) template const xcomplex<utype>& xcomplex<utype>::operator=(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ME20200107 BEGIN
  template <class utype> bool identical(const xcomplex<utype>& a, const xcomplex<utype>& b, utype _tol_) {
    return (identical(a.re, b.re, _tol_) && identical(a.im, b.im, _tol_));
  }
#define AST_TEMPLATE(utype) template bool identical(const xcomplex<utype>&, const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  template <class utype> bool isdifferent(const xcomplex<utype>& a, const xcomplex<utype>& b, utype _tol_) {
    return !identical(a, b, _tol_);
  }
#define AST_TEMPLATE(utype) template bool isdifferent(const xcomplex<utype>&, const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  template <class utype> bool isdifferent(const xcomplex<utype>& a, const xcomplex<utype>& b) {
    return !identical(a, b, (utype) _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
  }
#define AST_TEMPLATE(utype) template bool isdifferent(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  template <class utype> bool isequal(const xcomplex<utype>& a, const xcomplex<utype>& b, utype _tol_) {
    return identical(a, b, _tol_);
  }
#define AST_TEMPLATE(utype) template bool isequal(const xcomplex<utype>&, const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  template <class utype> bool isequal(const xcomplex<utype>& a, const xcomplex<utype>& b) {
    return identical(a, b, (utype) _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
  }
#define AST_TEMPLATE(utype) template bool isequal(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE
  // ME20200107 END

  // namespace aurostd
  //   template<class utype>                            // operator <<  xcomplex<>
  //   std::ostream& operator<< (std::ostream& buf,const xcomplex<utype>& x) {
  //     buf << "[" << x.re << "," << x.im << "]";
  //     return buf;
  //   }

  // ---------------------------------------------- operator ostream << xcomplex
  template <class utype> std::ostream& operator<<(std::ostream& os, const xcomplex<utype>& x) {
    return os << '(' << real(x) << ',' << imag(x) << ')';
  }
#define AST_TEMPLATE(utype) template std::ostream& operator<<(std::ostream&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator +=
  template <class utype> xcomplex<utype>& xcomplex<utype>::operator+=(const xcomplex<utype>& r) {
    this->re += r.re;
    this->im += r.im;
    return *this;
  }
#define AST_TEMPLATE(utype) template xcomplex<utype>& xcomplex<utype>::operator+=(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator -=
  template <class utype> xcomplex<utype>& xcomplex<utype>::operator-=(const xcomplex<utype>& r) {
    this->re -= r.re;
    this->im -= r.im;
    return *this;
  }
#define AST_TEMPLATE(utype) template xcomplex<utype>& xcomplex<utype>::operator-=(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator *=
  template <class utype> xcomplex<utype>& xcomplex<utype>::operator*=(const xcomplex<utype>& r) {
    utype freal = this->re * r.re - this->im * r.im;
    utype fimag = this->re * r.im + this->im * r.re;
    this->im = fimag;
    this->re = freal;
    return *this;
  }
#define AST_TEMPLATE(utype) template xcomplex<utype>& xcomplex<utype>::operator*=(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator *=
  template <class utype> xcomplex<utype>& xcomplex<utype>::operator/=(const xcomplex<utype>& r) {
    utype ar = abs(r.re);
    utype ai = abs(r.im);
    utype xr;
    utype ni;
    utype t;
    utype d;
    if (ar <= ai) {
      t = r.re / r.im;
      d = r.im * (1 + t * t);
      xr = (this->re * t + this->im) / d;
      ni = (this->im * t - this->re) / d;
    } else {
      t = r.im / r.re;
      d = r.re * (1 + t * t);
      xr = (this->re + this->im * t) / d;
      ni = (this->im - this->re * t) / d;
    }
    this->re = xr;
    this->im = ni;
    return *this;
  }
#define AST_TEMPLATE(utype) template xcomplex<utype>& xcomplex<utype>::operator/=(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator ==
  template <class utype> bool xcomplex<utype>::operator==(const xcomplex<utype>& r) {
    utype a = this->re * this->re + this->im + this->im;
    utype b = r.re * r.re + r.im * r.im;
    return a == b;
  }
#define AST_TEMPLATE(utype) template bool xcomplex<utype>::operator==(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator <
  template <class utype> bool xcomplex<utype>::operator<(const xcomplex<utype>& r) {
    utype a = this->re * this->re + this->im + this->im;
    utype b = r.re * r.re + r.im * r.im;
    return a < b;
  }
#define AST_TEMPLATE(utype) template bool xcomplex<utype>::operator<(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator >
  template <class utype> bool xcomplex<utype>::operator>(const xcomplex<utype>& r) {
    utype a = this->re * this->re + this->im + this->im;
    utype b = r.re * r.re + r.im * r.im;
    return a > b;
  }
#define AST_TEMPLATE(utype) template bool xcomplex<utype>::operator>(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator <=
  template <class utype> bool xcomplex<utype>::operator<=(const xcomplex<utype>& r) {
    utype a = this->re * this->re + this->im + this->im;
    utype b = r.re * r.re + r.im * r.im;
    return a <= b;
  }
#define AST_TEMPLATE(utype) template bool xcomplex<utype>::operator<=(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------------- operator >=
  template <class utype> bool xcomplex<utype>::operator>=(const xcomplex<utype>& r) {
    utype a = this->re * this->re + this->im + this->im;
    utype b = r.re * r.re + r.im * r.im;
    return a >= b;
  }
#define AST_TEMPLATE(utype) template bool xcomplex<utype>::operator>=(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------ function imag
  template <class utype> // removed inline
  utype imag(const xcomplex<utype>& x) {
    return x.imag();
  }
#define AST_TEMPLATE(utype) template utype imag(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------ function real
  template <class utype> // removed inline
  utype real(const xcomplex<utype>& x) {
    return x.real();
  }
#define AST_TEMPLATE(utype) template utype real(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ----------------------------------------------- operator xcomplex + xcomplex
  template <class utype> // removed inline
  xcomplex<utype> operator+(const xcomplex<utype>& x, const xcomplex<utype>& y) {
    return xcomplex<utype>(real(x) + real(y), imag(x) + imag(y));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator+(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------- operator xcomplex + utype
  template <class utype> // removed inline
  xcomplex<utype> operator+(const xcomplex<utype>& x, utype y) {
    return xcomplex<utype>(real(x) + y, imag(x));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator+(const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------- operator utype + xcomplex
  template <class utype> // removed inline
  xcomplex<utype> operator+(utype x, const xcomplex<utype>& y) {
    return xcomplex<utype>(x + real(y), imag(y));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator+(utype, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ----------------------------------------------- operator xcomplex - xcomplex
  template <class utype> // removed inline
  xcomplex<utype> operator-(const xcomplex<utype>& x, const xcomplex<utype>& y) {
    return xcomplex<utype>(real(x) - real(y), imag(x) - imag(y));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator-(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------- operator xcomplex - utype
  template <class utype> // removed inline
  xcomplex<utype> operator-(const xcomplex<utype>& x, utype y) {
    return xcomplex<utype>(real(x) - y, imag(x));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator-(const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------- operator utype - xcomplex
  template <class utype> // removed inline
  xcomplex<utype> operator-(utype x, const xcomplex<utype>& y) {
    return xcomplex<utype>(x - real(y), -imag(y));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator-(utype, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ----------------------------------------------- operator xcomplex * xcomplex
  template <class utype> // removed inline
  xcomplex<utype> operator*(const xcomplex<utype>& x, const xcomplex<utype>& y) {
    return xcomplex<utype>(real(x) * real(y) - imag(x) * imag(y), real(x) * imag(y) + imag(x) * real(y));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator*(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------- operator xcomplex * utype
  template <class utype> // removed inline
  xcomplex<utype> operator*(const xcomplex<utype>& x, utype y) {
    return xcomplex<utype>(real(x) * y, imag(x) * y);
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator*(const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------- operator utype * xcomplex
  template <class utype> // removed inline
  xcomplex<utype> operator*(utype x, const xcomplex<utype>& y) {
    return xcomplex<utype>(x * real(y), x * imag(y));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator*(utype, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ----------------------------------------------- operator xcomplex / xcomplex
  template <class utype> xcomplex<utype> operator/(const xcomplex<utype>& x, const xcomplex<utype>& y) {
    utype ar = abs(real(y));
    utype ai = abs(imag(y));
    utype xrr;
    utype nir;
    utype xri;
    utype nii;
    utype t;
    utype d;
    utype xr = abs(real(x));
    utype xi = abs(imag(x));
    if (ar <= ai) {
      t = real(y) / imag(y);
      d = imag(y) * (1 + t * t);
      xrr = xr * t / d;
      nir = -xr / d;
      xri = xi * t / d;
      nii = -xi / d;
    } else {
      t = imag(y) / real(y);
      d = real(y) * (1 + t * t);
      xrr = xr / d;
      nir = -xr * t / d;
      xri = xi / d;
      nii = -xi * t / d;
    }
    return xcomplex<utype>(xrr + nii, nir + xri);
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator/(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------- operator utype / xcomplex
  template <class utype> xcomplex<utype> operator/(utype x, const xcomplex<utype>& y) {
    utype ar = abs(real(y));
    utype ai = abs(imag(y));
    utype xr;
    utype ni;
    utype t;
    utype d;
    if (ar <= ai) {
      t = real(y) / imag(y);
      d = imag(y) * (1 + t * t);
      xr = x * t / d;
      ni = -x / d;
    } else {
      t = imag(y) / real(y);
      d = real(y) * (1 + t * t);
      xr = x / d;
      ni = -x * t / d;
    }
    return xcomplex<utype>(xr, ni);
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator/(utype, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  //  // ------------------------------------------------- operator utype / xcomplex
  //  template<class utype> xcomplex<utype>
  //    operator /(const utype x,const xcomplex<utype>& y) {
  //    return xcomplex<utype>(real(x) / y,imag(x) / y); //  MISSING
  //  }

  // ------------------------------------------------- operator xcomplex / utype
  template <class utype> xcomplex<utype> operator/(const xcomplex<utype>& x, utype y) {
    return xcomplex<utype>(real(x) / y, imag(x) / y);
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator/(const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  template <class utype1, class utype2> xcomplex<utype1> operator/(const xcomplex<utype1>& x, utype2 y) {
    return xcomplex<utype1>(real(x) / (utype1) y, imag(x) / (utype1) y);
  }

  // -------------------------------------------------------- operator +xcomplex
  template <class utype> // removed inline
  xcomplex<utype> operator+(const xcomplex<utype>& x) {
    return x;
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator+(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // -------------------------------------------------------- operator -xcomplex
  template <class utype> // removed inline
  xcomplex<utype> operator-(const xcomplex<utype>& x) {
    return xcomplex<utype>(-real(x), -imag(x));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> operator-(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ---------------------------------------------- operator xcomplex == xcomplex
  template <class utype> // removed inline
  bool operator==(const xcomplex<utype>& x, const xcomplex<utype>& y) {
    return real(x) == real(y) && imag(x) == imag(y);
  }
#define AST_TEMPLATE(utype) template bool operator==(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------ operator xcomplex == utype
  template <class utype> // removed inline
  bool operator==(const xcomplex<utype>& x, utype y) {
    return real(x) == y && imag(x) == 0;
  }
#define AST_TEMPLATE(utype) template bool operator==(const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------ operator utype == xcomplex
  template <class utype> // removed inline
  bool operator==(utype x, const xcomplex<utype>& y) {
    return x == real(y) && imag(y) == 0;
  }
#define AST_TEMPLATE(utype) template bool operator==(utype, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ---------------------------------------------- operator xcomplex != xcomplex
  template <class utype> // removed inline
  bool operator!=(const xcomplex<utype>& x, const xcomplex<utype>& y) {
    return real(x) != real(y) || imag(x) != imag(y);
  }
#define AST_TEMPLATE(utype) template bool operator!=(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------ operator xcomplex != utype
  template <class utype> // removed inline
  bool operator!=(const xcomplex<utype>& x, utype y) {
    return real(x) != y || imag(x) != 0;
  }
#define AST_TEMPLATE(utype) template bool operator!=(const xcomplex<utype>&, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------ operator utype != xcomplex
  template <class utype> // removed inline
  bool operator!=(utype x, const xcomplex<utype>& y) {
    return x != real(y) || imag(y) != 0;
  }
#define AST_TEMPLATE(utype) template bool operator!=(utype, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------- function abs
  extern "C" double hypot(double, double) __xprototype;
  template <class utype> // removed inline
  utype abs(const xcomplex<utype>& x) {
    return hypot(real(x), imag(x));
  }
#define AST_TEMPLATE(utype) template utype abs(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------- function arg
  template <class utype> // removed inline
  utype arg(const xcomplex<utype>& x) {
    return atan2(imag(x), real(x));
  }
#define AST_TEMPLATE(utype) template utype arg(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ----------------------------------------------------------- function polar
  template <class utype> // removed inline
  xcomplex<utype> polar(utype r, utype t) {
    return xcomplex<utype>(r * std::cos(t), r * std::sin(t));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> polar(utype, utype);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------ function conj
  template <class utype> // removed inline
  xcomplex<utype> conj(const xcomplex<utype>& x) {
    return xcomplex<utype>(real(x), -imag(x));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> conj(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------ function conj
  // ME20180904 -- needed to compile conjugate transpose in xmatrix
  template <class utype> utype conj(const utype& x) {
    return x;
  }
#define AST_TEMPLATE(utype) template utype conj(const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // ------------------------------------------------------------ function norm
  template <class utype> // removed inline
  utype norm(const xcomplex<utype>& x) {
    return real(x) * real(x) + imag(x) * imag(x);
  }
#define AST_TEMPLATE(utype) template utype norm(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------- function cos
  template <class utype> xcomplex<utype> cos(const xcomplex<utype>& x) {
    return xcomplex<utype>(std::cos(real(x)) * std::cosh(imag(x)), -std::sin(real(x)) * std::sinh(imag(x)));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> cos(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------ function cosh
  template <class utype> xcomplex<utype> cosh(const xcomplex<utype>& x) {
    return xcomplex<utype>(std::cosh(real(x)) * std::cos(imag(x)), std::sinh(real(x)) * std::sin(imag(x)));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> cosh(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------- function exp
  template <class utype> xcomplex<utype> exp(const xcomplex<utype>& x) {
    return polar(utype(std::exp(real(x))), imag(x));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> exp(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------- function log
  template <class utype> xcomplex<utype> log(const xcomplex<utype>& x) {
    return xcomplex<utype>(std::log(abs(x)), arg(x));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> log(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // --------------------------------------------- function xcomplex pow xcomplex
  template <class utype> xcomplex<utype> pow(const xcomplex<utype>& x, const xcomplex<utype>& y) {
    utype logr = std::log(abs(x));
    utype t = arg(x);
    return polar(utype(std::exp(logr * real(y) - imag(y) * t)), utype(imag(y) * logr + real(y) * t));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> pow(const xcomplex<utype>&, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------- function xcomplex pow int
  template <class utype> xcomplex<utype> pow(const xcomplex<utype>& xin, int y) {
    if (y == 0) {
      return xcomplex<utype>(1.0);
    }
    xcomplex<utype> r(1.0);
    xcomplex<utype> x(xin);
    if (y < 0) {
      y = -y;
      x = ((utype) 1.0) / x;
    }
    for (;;) {
      if (y & 1) {
        r *= x;
      }
      if (y >>= 1) {
        x *= x;
      } else {
        return r;
      }
    }
  }
// #define AST_TEMPLATE(utype) template xcomplex<utype> pow(const xcomplex<utype>&,int);
//     AST_GEN_1(AST_UTYPE_FLOAT)
// #undef AST_TEMPLATE

  // ----------------------------------------------- function xcomplex pow utype
  template <class utype> xcomplex<utype> pow(const xcomplex<utype>& x, utype y) {
    return exp(utype(y) * log(x));
  }
// #define AST_TEMPLATE(utype) template xcomplex<utype> pow(const xcomplex<utype>&,utype);
//     AST_GEN_1(AST_UTYPE_FLOAT)
// #undef AST_TEMPLATE

  // ----------------------------------------------- function utype pow xcomplex
  template <class utype> xcomplex<utype> pow(utype x, const xcomplex<utype>& y) {
    return exp(y * utype(std::log(x)));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> pow(utype, const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------- function sin
  template <class utype> xcomplex<utype> sin(const xcomplex<utype>& x) {
    return xcomplex<utype>(std::sin(real(x)) * std::cosh(imag(x)), std::cos(real(x)) * std::sinh(imag(x)));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> sin(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------ function sinh
  template <class utype> xcomplex<utype> sinh(const xcomplex<utype>& x) {
    return xcomplex<utype>(std::sinh(real(x)) * std::cos(imag(x)), std::cosh(real(x)) * std::sin(imag(x)));
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> sinh(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ------------------------------------------------------------ function sqrt
  template <class utype> xcomplex<utype> sqrt(const xcomplex<utype>& x) {
    utype r = abs(x);
    utype xr;
    utype ni;
    if (r == 0.0) {
      xr = ni = r;
    } else if (real(x) > 0) {
      xr = std::sqrt(0.5 * (r + real(x)));
      ni = imag(x) / xr / 2;
    } else {
      ni = std::sqrt(0.5 * (r - real(x)));
      if (imag(x) < 0) {
        ni = -ni;
      }
      xr = imag(x) / ni / 2;
    }
    return xcomplex<utype>(xr, ni);
  }
#define AST_TEMPLATE(utype) template xcomplex<utype> sqrt(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ME20180907
  //  ------------------------------------------------------------ function magnitude square
  template <class utype> utype magsqr(const xcomplex<utype>& x) {
    return (x.re * x.re + x.im * x.im);
  }
#define AST_TEMPLATE(utype) template utype magsqr(const xcomplex<utype>&);
  AST_GEN_1(AST_UTYPE_FLOAT)
#undef AST_TEMPLATE

  // ---------------------------------------------- operator istream >> xcomplex
  // WARNING SC FIGURE OUT WHERE ipfx0() and isfx() went !

  // template<class utype> istream&
  // operator >>(istream& is,xcomplex<utype>& x) {
  // utype re,im=0;
  // char ch=0;
  // if(is.ipfx0()) {
  // if(is.peek() == '(')
  // is >> ch;
  // is >> re;
  // if(ch == '(') {
  // is >> ch;
  // if(ch == ',')
  // is >> im >> ch;
  // }
  // }
  // is.isfx();
  // if(ch != 0 && ch != ')')
  // is.setstate(std::ios_base::failbit);
  // else if(is.good())
  // x=xcomplex<utype>(re,im);
  // return is;
  // }
  //  ---------------------------------------------- operator ostream << xcomplex

  // ----------------------------------------------------------------------------
  // ----------------------------------------- implementation for extra data type

} // namespace aurostd

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
