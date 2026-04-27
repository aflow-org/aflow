// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XSCALAR_H_
#define _AUROSTD_XSCALAR_H_

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <sys/types.h>

#include "aurostd_defs.h"
#include "aurostd_xcomplex.h"

#define _AUROSTD_XSCALAR_DEFAULT_SIZE_ 3
#define _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ 1.0e-6
#define _AUROSTD_XSCALAR_TOLERANCE_ROUNDOFF_ 1.0e-6

#define _AUROSTD_XSCALAR_TOLERANCE_INTEGER_ 1.0e-2 // DX20201217

namespace aurostd {
  // namespace aurostd
  bool _ishex(const std::string&) __xprototype;
  template <class utype> bool _isodd(utype) __xprototype;
  template <class utype> bool _iseven(utype) __xprototype;
  template <class utype> bool _isreal(utype) __xprototype;
  template <class utype> bool _iscomplex(utype) __xprototype;
  template <class utype> size_t _size(utype) __xprototype;
  template <class utype> size_t _size(xcomplex<utype>) __xprototype;
  template <class utype> utype abs(utype) __xprototype;
  template <class utype> double sqrt(utype) __xprototype;
  long double sqrt(long double) __xprototype;
  template <class utype> double exp(utype) __xprototype;
  long double exp(long double) __xprototype;
  template <class utype> utype sign(utype) __xprototype;
  template <class utype> utype mod(utype, utype) __xprototype;
  template <class utype> utype mod_floored(utype, utype) __xprototype;
  template <class utype> utype nint(utype) __xprototype;
  template <class utype> void _GCD(int a, int b, int& gcd, int& x, int& y); // CO20180409  //CO20191112 - extended GCD, get Bezout coefficients
  template <class utype> void _GCD(int a, int b, int& gcd); // CO20180409
  void GCD(int a, int b, int& gcd, int& x, int& y); // CO20191201
  void GCD(int a, int b, int& gcd); // CO20191201
  void GCD(uint a, uint b, uint& gcd, uint& x, uint& y);  // CO20191201
  void GCD(uint a, uint b, uint& gcd);  // CO20191201
  void GCD(long int a, long int b, long int& gcd, long int& x, long int& y);  // CO20191201
  void GCD(long int a, long int b, long int& gcd);  // CO20191201
  void GCD(unsigned long int a, unsigned long int b, unsigned long int& gcd, unsigned long int& x, unsigned long int& y); // CO20191201
  void GCD(unsigned long int a, unsigned long int b, unsigned long int& gcd); // CO20191201
  void GCD(long long int a, long long int b, long long int& gcd, long long int& x, long long int& y); // CO20191201
  void GCD(long long int a, long long int b, long long int& gcd); // CO20191201
  void GCD(unsigned long long int a, unsigned long long int b, unsigned long long int& gcd, unsigned long long int& x, unsigned long long int& y);  // CO20191201
  void GCD(unsigned long long int a, unsigned long long int b, unsigned long long int& gcd);  // CO20191201
  void GCD(float a, float b, float& gcd, float& x, float& y, float tolerance = 0.01);  // CO20191201
  void GCD(float a, float b, float& gcd, float tolerance = 0.01);  // CO20191201
  void GCD(double a, double b, double& gcd, double& x, double& y, double tolerance = 0.01);  // CO20191201
  void GCD(double a, double b, double& gcd, double tolerance = 0.01);  // CO20191201
  void GCD(long double a, long double b, long double& gcd, long double& x, long double& y, long double tolerance = 0.01);  // CO20191201
  void GCD(long double a, long double b, long double& gcd, long double tolerance = 0.01);  // CO20191201
  int LCM(int a, int b); // CO20190520

  template <class utype> bool isinteger(utype x) __xprototype;
  template <class atype, class btype> bool isinteger(atype x, btype tolerance) __xprototype;  // CO20191201
  template <class utype> bool iszero(utype a) __xprototype;
  template <class atype, class btype> bool iszero(atype a, btype tol) __xprototype;  // CO20191201

  std::string dbl2frac(double a, bool sign_prefix = true); // DX20190724
  double frac2dbl(const std::string& str); // DX20200313
  void double2fraction(const double& input_double, int& numerator, int& denominator, double tol_diff = AUROSTD_IDENTITY_TOL, double tol_remainder = 1e-2); // DX20210908
  int getNumeratorContinuedFractions(int& p, const int& n, std::vector<int>& fraction_sequence); // DX20210908
  int getDenominatorContinuedFractions(int& q, const int& n, std::vector<int>& fraction_sequence); // DX20210908

  template <class utype> utype fact(utype) __xprototype;
  template <class utype> utype factorial(utype) __xprototype;
  template <class utype> utype angle(utype, utype, utype, utype, utype, utype) __xprototype;
  template <class utype> utype angle(utype, utype, utype, utype) __xprototype;
  template <class utype> utype modulus(utype, utype, utype) __xprototype;
  template <class utype> utype modulus(utype, utype) __xprototype;

  double ln(double);
  float lnf(float);
  //  long double lnl(long double);
  float ln(float);
  //  long double ln(long double);
  double log(double);
  float logf(float);
  float log(float);
  //  long double logl(long double);
  //  long double log(long double);
  double log10(double);
  float log10f(float);
  float log10(float);
  // long double log10l(long double);
  // long double log10(long double);
} // namespace aurostd

// ----------------------------------------------------------------------------
// _isfloat  _isfloat  _isfloat  _isfloat  _isfloat
namespace aurostd {
  // namespace aurostd
  bool _isfloat(bool) __xprototype;
  bool _isfloat(char) __xprototype;
  bool _isfloat(uint) __xprototype;
  bool _isfloat(int) __xprototype;
  bool _isfloat(long int) __xprototype;
  bool _isfloat(unsigned long int) __xprototype;
  bool _isfloat(long long int) __xprototype;
  bool _isfloat(unsigned long long int) __xprototype;
  bool _isfloat(float) __xprototype;
  bool _isfloat(double) __xprototype;
  bool _isfloat(long double) __xprototype;
  bool _isfloat(xcomplex<float>) __xprototype;
  bool _isfloat(xcomplex<double>) __xprototype;
  bool _isfloat(xcomplex<long double>) __xprototype;
  bool isfloat(const std::string& in); // CO20180729
} // namespace aurostd

// ----------------------------------------------------------------------------
// _iscomplex  _iscomplex  _iscomplex  _iscomplex  _iscomplex
namespace aurostd {
  // namespace aurostd
  bool _iscomplex(bool) __xprototype;
  bool _iscomplex(char) __xprototype;
  bool _iscomplex(uint) __xprototype;
  bool _iscomplex(int) __xprototype;
  bool _iscomplex(long int) __xprototype;
  bool _iscomplex(unsigned long int) __xprototype;  // CO20191201
  bool _iscomplex(long long int) __xprototype;
  bool _iscomplex(unsigned long long int) __xprototype;
  bool _iscomplex(float) __xprototype;
  bool _iscomplex(double) __xprototype;
  bool _iscomplex(long double) __xprototype;
  bool _iscomplex(xcomplex<float>) __xprototype;
  bool _iscomplex(xcomplex<double>) __xprototype;
  bool _iscomplex(xcomplex<long double>) __xprototype;
} // namespace aurostd

// ----------------------------------------------------------------------------
// _isreal  _isreal  _isreal  _isreal  _isreal
namespace aurostd {
  // namespace aurostd
  bool _isreal(bool) __xprototype;
  bool _isreal(char) __xprototype;
  bool _isreal(uint) __xprototype;
  bool _isreal(int) __xprototype;
  bool _isreal(long int) __xprototype;
  bool _isreal(unsigned long int) __xprototype; // CO20191201
  bool _isreal(long long int) __xprototype;
  bool _isreal(unsigned long long int) __xprototype;
  bool _isreal(float) __xprototype;
  bool _isreal(double) __xprototype;
  bool _isreal(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _isreal(xcomplex<float>) __xprototype;
  bool _isreal(xcomplex<double>) __xprototype;
  bool _isreal(xcomplex<long double>) __xprototype;
#endif
} // namespace aurostd

// ----------------------------------------------------------------------------
// _isreal  _isreal  _isreal  _isreal  _isreal
namespace aurostd {
  // namespace aurostd
  bool _isreal(bool) __xprototype;
  bool _isreal(char) __xprototype;
  bool _isreal(uint) __xprototype;
  bool _isreal(int) __xprototype;
  bool _isreal(long int) __xprototype;
  bool _isreal(unsigned long int) __xprototype; // CO20191201
  bool _isreal(long long int) __xprototype;
  bool _isreal(unsigned long long int) __xprototype;
  bool _isreal(float) __xprototype;
  bool _isreal(double) __xprototype;
  bool _isreal(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _isreal(xcomplex<float>) __xprototype;
  bool _isreal(xcomplex<double>) __xprototype;
  bool _isreal(xcomplex<long double>) __xprototype;
#endif
} // namespace aurostd

// ----------------------------------------------------------------------------
// _iseven  _iseven  _iseven  _iseven  _iseven
namespace aurostd {
  // namespace aurostd
  bool _iseven(bool) __xprototype;
  bool _iseven(char) __xprototype;
  bool _iseven(int) __xprototype;
  bool _iseven(uint) __xprototype;
  bool _iseven(float) __xprototype;
  bool _iseven(double) __xprototype;
  bool _iseven(long int) __xprototype;
  bool _iseven(long long int) __xprototype;
  bool _iseven(unsigned long long int) __xprototype;
  bool _iseven(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _iseven(xcomplex<float>) __xprototype;
  bool _iseven(xcomplex<double>) __xprototype;
  bool _iseven(xcomplex<long double>) __xprototype;
#endif
} // namespace aurostd

// ----------------------------------------------------------------------------
// _isodd  _isodd  _isodd  _isodd  _isodd
namespace aurostd {
  // namespace aurostd
  bool _isodd(bool) __xprototype;
  bool _isodd(char) __xprototype;
  bool _isodd(int) __xprototype;
  bool _isodd(uint) __xprototype;
  bool _isodd(float) __xprototype;
  bool _isodd(double) __xprototype;
  bool _isodd(long int) __xprototype;
  bool _isodd(long long int) __xprototype;
  bool _isodd(unsigned long long int) __xprototype;
  bool _isodd(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _isodd(xcomplex<float>) __xprototype;
  bool _isodd(xcomplex<double>) __xprototype;
  bool _isodd(xcomplex<long double>) __xprototype;
#endif
} // namespace aurostd

// ----------------------------------------------------------------------------
// _real _real _real _real _real _real _real _real
namespace aurostd {
  // namespace aurostd
  bool _real(bool) __xprototype;
  char _real(char) __xprototype;
  uint _real(uint) __xprototype;
  int _real(int) __xprototype;
  long int _real(long int) __xprototype;
  unsigned long int _real(unsigned long int) __xprototype;  // CO20191201
  long long int _real(long long int) __xprototype;
  unsigned long long int _real(unsigned long long int) __xprototype;
  float _real(float) __xprototype;
  double _real(double) __xprototype;
  long double _real(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  float _real(xcomplex<float>) __xprototype;
  double _real(xcomplex<double>) __xprototype;
  long double _real(xcomplex<long double>) __xprototype;
#endif
} // namespace aurostd

// ----------------------------------------------------------------------------
// signnozero signnozero signnozero signnozero signnozero signnozero signnozero signnozero
namespace aurostd {
  // namespace aurostd
  bool signnozero(bool) __xprototype;
  char signnozero(char) __xprototype;
  int signnozero(int) __xprototype;
  uint signnozero(uint) __xprototype;
  float signnozero(float) __xprototype;
  double signnozero(double) __xprototype;
  long int signnozero(long int) __xprototype;
  long long int signnozero(long long int) __xprototype;
  long double signnozero(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  xcomplex<float> signnozero(xcomplex<float>) __xprototype;
  xcomplex<double> signnozero(xcomplex<double>) __xprototype;
  xcomplex<long double> signnozero(xcomplex<long double>) __xprototype;
#endif
} // namespace aurostd

// ----------------------------------------------------------------------------
// nint nint nint nint nint nint nint nint
namespace aurostd {
  // namespace aurostd
  bool nint(bool) __xprototype;
  char nint(char) __xprototype;
  int nint(int) __xprototype;
  uint nint(uint) __xprototype;
  uint nint(uint) __xprototype;
  float nint(float) __xprototype;
  double nint(double) __xprototype;
  long int nint(long int) __xprototype;
  long long int nint(long long int) __xprototype;
  long long unsigned int nint(long long unsigned int) __xprototype;
  long double nint(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  xcomplex<float> nint(xcomplex<float>) __xprototype;
  xcomplex<double> nint(xcomplex<double>) __xprototype;
  xcomplex<long double> nint(xcomplex<long double>) __xprototype;
#endif
} // namespace aurostd

// ----------------------------------------------------------------------------
// abs  abs  abs  abs  abs
namespace aurostd {
  // namespace aurostd
  // ABS(X)
  // int abs(int x) __xprototype;
  char abs(char x) __xprototype;
  int abs(int x) __xprototype;
  uint abs(uint x) __xprototype;
  float abs(float x) __xprototype;
  double abs(double x) __xprototype;
  long int abs(long int x) __xprototype;
  long long int abs(long long int x) __xprototype;
  unsigned long int abs(unsigned long int x) __xprototype;
  unsigned long long int abs(unsigned long long int x) __xprototype;
  long double abs(long double x) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  float abs(xcomplex<float> x) __xprototype;
  double abs(xcomplex<double> x) __xprototype;
  long double abs(xcomplex<long double> x) __xprototype;
#endif
} // namespace aurostd

//----------------------------------------------------------------------------
// pow  pow  pow  pow  pow
namespace aurostd {
  template <class utype> double pow(utype, utype);
  long double pow(long double, long double);
}

//----------------------------------------------------------------------------
// sin  sin  sin  sin  sin
namespace aurostd {
  template <class utype> double sin(utype);
    long double sin(long double);
} 

//----------------------------------------------------------------------------
// cos  cos  cos  cos  cos
namespace aurostd {
  template <class utype> double cos(utype);
  long double cos(long double);
} 

//----------------------------------------------------------------------------
// atan2  atan2  atan2  atan2  atan2
namespace aurostd {
  template <class utype> double atan2(utype, utype);
  long double atan2(long double, long double);
}

//--------------------------------------------------------------- round
namespace aurostd {
  double round(double x, uint digits = 0); // SD20220603
  int roundDouble(double doub, int multiple, bool up);  // CO20220624
  bool greaterEqualZero(double val);  // CO20220624
  bool lessEqualZero(double val); // CO20220624
  bool notPositive(double val, bool soft_cutoff, double tol);  // CO20220624
  bool notNegative(double val, bool soft_cutoff, double tol);  // CO20220624
  bool zeroWithinTol(double val, double tol); // CO20220624
  bool nonZeroWithinTol(double val, double tol);  // CO20220624
} // namespace aurostd

//--------------------------------------------------------------- isequal
namespace aurostd {
  // with const utype&
  template <class utype> bool identical(const utype&, const utype&, const utype&) __xprototype;
  template <class utype> bool identical(const utype&, const utype&) __xprototype;
  bool identical(const bool a, const bool b); // SD20220705
  bool identical(const char a, const char b); // SD20220705
  bool identical(const std::string& a, const std::string& b); // SD20220705
  template <class utype> bool isdifferent(const utype&, const utype&, const utype&) __xprototype;
  template <class utype> bool isdifferent(const utype&, const utype&) __xprototype;
  bool isdifferent(const bool a, const bool b); // SD20220705
  bool isdifferent(const char a, const char b); // SD20220705
  bool isdifferent(const std::string& a, const std::string& b); // SD20220705
  template <class utype> bool isequal(const utype&, const utype&, const utype&) __xprototype;
  template <class utype> bool isequal(const utype&, const utype&) __xprototype;
  bool isequal(const bool a, const bool b); // SD20220705
  bool isequal(const char a, const char b); // SD20220705
  bool isequal(const std::string& a, const std::string& b); // SD20220705
  // with utype
  // template<class utype> bool identical(utype,utype,utype) __xprototype;
  // template<class utype> bool identical(utype,utype) __xprototype;
  // template<class utype> bool isdifferent(utype,utype,utype) __xprototype;
  // template<class utype> bool isdifferent(utype,utype) __xprototype;
  // template<class utype> bool isequal(utype,utype,utype) __xprototype;
  // template<class utype> bool isequal(utype,utype) __xprototype;
} // namespace aurostd

//--------------------------------------------------------------- extra min/max // __XEXTRA_MINMAX_CPP
namespace aurostd {
  // namespace aurostd
  template <class utype> utype min(utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype min(utype, utype, utype, utype, utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template <class utype> utype max(utype, utype, utype, utype, utype, utype, utype, utype, utype, utype, utype, utype) __xprototype;  // __XEXTRA_MINMAX_CPP
} // namespace aurostd

namespace aurostd {
  template <class utype> utype _roundoff(const utype& x, utype tolerance);
  int roundoff(int x, int tolerance);
  long roundoff(long x, long tolerance);
  uint roundoff(uint x, uint tolerance);
  float roundoff(float x, float tolerance);
  double roundoff(double x, double tolerance);
  long long int roundoff(long long int x, long long int tolerance);
  unsigned long int roundoff(unsigned long int x, unsigned long int tolerance);
  unsigned long long int roundoff(unsigned long long int x, unsigned long long int tolerance);
  long double roundoff(long double x, long double tolerance);
} // namespace aurostd

namespace aurostd {
  int boundary_conditions_periodic(int lrows, int urows, int i);  // CO20190419 - taken from xvector BOUNDARY_CONDITIONS_PERIODIC
}

namespace aurostd {
  uint powint(uint x, uint exp); // CO20191201
  int powint(int x, uint exp); // CO20191201
  long int powint(long int x, uint exp); // CO20191201
  unsigned long int powint(unsigned long int x, uint exp); // CO20191201
  long long int powint(long long int x, uint exp); // CO20191201
  unsigned long long int powint(unsigned long long int x, uint exp); // CO20191201
} // namespace aurostd

// AS20200513 BEGIN
namespace aurostd {
  double FermiDirac(double E, double mu, double T);
}
// AS20200513 END

// CO20201111 BEGIN
namespace aurostd {
  template <class utype> utype nCk(utype n, utype k);
}
// CO20201111 END

// CO20201111 - BEGIN
namespace aurostd {
  bool isNaN(double d);
}
// CO20201111 - END

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XSCALAR_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
