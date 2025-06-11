// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XSCALAR_H_
#define _AUROSTD_XSCALAR_H_
#define _AUROSTD_XSCALAR_DEFAULT_SIZE_ 3

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
#ifndef __xprototype
#define __xprototype __attribute__((const))
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

// ----------------------------------------------------------------------------
// _AUROSTD_XSCALAR_PROTOTYPES_

// ----------------------------------------------------------------------------
// ---------------------------------------------------------- XXXXXXXXXXXXXXXXX

// ----------------------------------------------------------------------------
// ---------------------------------------------------------- physics constants

#define rad2deg (180.0/3.14159265358979323846)
#define deg2rad (3.14159265358979323846/180.0)
#define angstrom2bohr (1/0.529177249)
#define bohr2angstrom (0.529177249)

//ME20181020 - check out https://physics.nist.gov/cuu/Constants/index.html
//#define PI                    3.14159265359
#define PI                      3.14159265358979323846
#define pi                      PI
#define TWOPI                   6.28318530717958647692
#define EULERSNUMBER            2.71828182845904523536
#define RTPI                    (sqrt(PI))
#define C_VACUUM                2.99792458E+8                   // m/s
#define EPS_VACUUM              8.854187817E-12                 // C/(N-m2)
#define MU_VACUUM               (4.0*PI*1.0E-7)                 // N/A^2 (T^2m^3/J)
#define AMU2KILOGRAM            1.66054E-27
#define KILOGRAM2AMU            6.022137E+26
#define E_ELECTRON              1.60217662E-19                  // C
#define eV2J                    E_ELECTRON                      // 1eV=E_ELECTRON J //CO20201111
#define J2eV                    (1.0/E_ELECTRON)                //CO20201111
#define PLANCKSCONSTANT_h       6.62607E-34                     // J*s
#define PLANCKSCONSTANT_hbar    1.0545718E-34                   // J*s
#define PLANCKSCONSTANTEV_h     (PLANCKSCONSTANT_h/E_ELECTRON)  // eV*s
#define PLANCKSCONSTANTEV_hbar  (PLANCKSCONSTANTEV_h/TWOPI)     // eV*s
#define KBOLTZ                  1.3806504E-23                   // J/K
#define KBOLTZEV                (KBOLTZ/E_ELECTRON)             // eV/K //8.617343E-5
#define eV2K                    (11604.505)                     // 1eV=11604.505 K
#define meV2K                   (11.604505)                     // 1meV=11.604505 K
#define mol2atom                6.0221408E23                    // 1mol=6.022e23 atoms    //CO20180329
#define atom2mol                (1.0/6.0221408E23)              //CO20201111
#define eVatom2kJmol            (E_ELECTRON*mol2atom/1.0e3)     // 1eV/atom=96.5kJ/mol    //CO20180329
#define meVatom2kJmol           (eVatom2kJmol/1.0e3)            // 1meV/atom=0.0965kJ/mol //CO20180329
#define hartree2eV              27.2113862459                   // 1hartree=27.211eV      //ME20200206
#define kcal2eV                 4.336443203200000E-002          // 1(kcal/mol) = 4.33644E-2 eV
#define eV2kcal                 (1.0/kcal2eV)

//ME20200107 - (A)APL conversion factors
#define THz2Hz                        1E12
#define Hz2THz                        (1.0/THz2Hz)
#define au2Hz                         aurostd::sqrt(E_ELECTRON*1E20/AMU2KILOGRAM)  // eV/(A^2 amu) -> Hz
#define au2rcm                        (au2Hz/(100*C_VACUUM))                         // eV/(A^2 amu) -> cm^-1
#define au2eV                         (au2Hz*PLANCKSCONSTANTEV_h)                    // eV/(A^2 amu) -> eV
#define eV2Hz                         (1.0/PLANCKSCONSTANTEV_h)
#define eV2rcm                        (1.0/(PLANCKSCONSTANTEV_h*100*C_VACUUM))
#define au2nmTHz                      ((E_ELECTRON*Hz2THz*Hz2THz*1E18)/(0.1 * AMU2KILOGRAM))  // eV/(A amu) -> nm * THz^2
#define PLANCKSCONSTANT_h_THz         (PLANCKSCONSTANT_h*THz2Hz) // J/THz
#define PLANCKSCONSTANT_hbar_THz      (PLANCKSCONSTANT_hbar*THz2Hz) // J/THz
#define PLANCKSCONSTANTAMU_hbar_THz   (PLANCKSCONSTANTEV_hbar*THz2Hz*(10*au2nmTHz))  // amu A^2 THz
#define BEfactor_hbar_THz             (PLANCKSCONSTANTEV_hbar/(KBOLTZEV*Hz2THz))  // hbar/kB in K/THz
#define BEfactor_h_THz                (PLANCKSCONSTANTEV_h/(KBOLTZEV*Hz2THz))  // h/kB in K/THz

//AS20200427 - QHA-related conversion factors
#define eV2GPa (E_ELECTRON*1e21)    // [eV/A^3] --> [GPa]
#define GPa2eV (1.0/eV2GPa)         // [GPa] --> [eV/A^3]
#define eV2kBar (eV2GPa*10)         // [eV/A^3] --> [kBar]
#define kBar2eV (1.0/eV2kBar)       // [kBar] --> [eV/A^3]
#define atm2Pa 101325

//DX20210111 - GFA factors
#define TEMPERATURE_ROOM 300.0               // K
#define kBT_ROOM (KBOLTZEV*TEMPERATURE_ROOM) // 0.025

// ----------------------------------------------------------------------------
// ------------------------------------------------------------------ constants

#define _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ 1.0e-6
#define _AUROSTD_XSCALAR_TOLERANCE_ROUNDOFF_ 1.0e-6

#define _AUROSTD_XSCALAR_TOLERANCE_INTEGER_ 1.0e-2 //DX20201217

// ----------------------------------------------------------------------------
// ------------------------- primitives for template<class utype> xscalar<utype>

// ----------------------------------------------------------------------------
// ------------------------------------------------- unary and binary operators

//--------------------------------------------------------- extra on data types

// namespace aurostd {
//   double atof(string str);
//   int atoi(string str);
//   long atol(string str);
//   // long long atoll(string str);
//   // long long atoq(string str);
// }

namespace aurostd {
  // namespace aurostd
  //[CO20180729 OBSOLETE]template<class utype> bool _isfloat(utype) __xprototype;
  bool _ishex(const string&) __xprototype;
  template<class utype> bool _isodd(utype) __xprototype;
  template<class utype> bool _iseven(utype) __xprototype;
  template<class utype> bool _isreal(utype) __xprototype;
  template<class utype> bool _iscomplex(utype) __xprototype;
  template<class utype> int _size(utype) __xprototype;
  template<class utype> utype abs(utype) __xprototype;
  template<class utype> utype sqrt(utype) __xprototype;
  template<class utype> utype exp(utype) __xprototype; //SD20220512
  template<class utype> utype sign(utype) __xprototype;
  template<class utype> utype mod(utype,utype) __xprototype;
  template<class utype> utype mod_floored(utype,utype) __xprototype;
  template<class utype> utype nint(utype) __xprototype;
  template<class utype> void _GCD(int a,int b, int& gcd, int& x, int& y); //CO20180409  //CO20191112 - extended GCD, get Bezout coefficients
  template<class utype> void _GCD(int a,int b, int& gcd); //CO20180409
  void GCD(int a,int b,int& gcd,int& x,int& y); //CO20191201
  void GCD(int a,int b,int& gcd); //CO20191201
  void GCD(uint a,uint b,uint& gcd,uint& x,uint& y);  //CO20191201
  void GCD(uint a,uint b,uint& gcd);  //CO20191201
  void GCD(long int a,long int b,long int& gcd,long int& x,long int& y);  //CO20191201
  void GCD(long int a,long int b,long int& gcd);  //CO20191201
  void GCD(unsigned long int a,unsigned long int b,unsigned long int& gcd,unsigned long int& x,unsigned long int& y); //CO20191201
  void GCD(unsigned long int a,unsigned long int b,unsigned long int& gcd); //CO20191201
  void GCD(long long int a,long long int b,long long int& gcd,long long int& x,long long int& y); //CO20191201
  void GCD(long long int a,long long int b,long long int& gcd); //CO20191201
  void GCD(unsigned long long int a,unsigned long long int b,unsigned long long int& gcd,unsigned long long int& x,unsigned long long int& y);  //CO20191201
  void GCD(unsigned long long int a,unsigned long long int b,unsigned long long int& gcd);  //CO20191201
  void GCD(float a,float b,float& gcd,float& x,float& y,float tolerance=0.01);  //CO20191201
  void GCD(float a,float b,float& gcd,float tolerance=0.01);  //CO20191201
  void GCD(double a,double b,double& gcd,double& x,double& y,double tolerance=0.01);  //CO20191201
  void GCD(double a,double b,double& gcd,double tolerance=0.01);  //CO20191201
  void GCD(long double a,long double b,long double& gcd,long double& x,long double& y,long double tolerance=0.01);  //CO20191201
  void GCD(long double a,long double b,long double& gcd,long double tolerance=0.01);  //CO20191201
  int LCM(int a,int b); //CO20190520

  template<class utype> bool _isinteger(utype,utype=(utype)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_) __xprototype;  //CO20191201
  //bool isinteger(bool x,bool tolerance=(bool)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_); //CO20191201
  //bool isinteger(char x,char tolerance=(char)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_); //CO20191201
  bool isinteger(uint x,uint tolerance=(uint)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_); //CO20191201  //CO20191201 - obvious but define for boot
  bool isinteger(int x,int tolerance=(int)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_);  //CO20191201  //CO20191201 - obvious but define for boot
  bool isinteger(long int x,long int tolerance=(long)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_); //CO20191201  //CO20191201 - obvious but define for boot
  bool isinteger(unsigned long int x,unsigned long int tolerance=(unsigned long)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_); //CO20191201  //CO20191201 - obvious but define for boot
  bool isinteger(long long int x,long long int tolerance=(long long int)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_);  //CO20191201  //CO20191201 - obvious but define for boot
  bool isinteger(unsigned long long int x,unsigned long long int tolerance=(unsigned long long int)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_); //CO20191201  //CO20191201 - obvious but define for boot
  bool isinteger(float x,float tolerance=(float)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_);  //CO20191201
  bool isinteger(double x,double tolerance=(double)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_); //CO20191201
  bool isinteger(long double x,long double tolerance=(long double)_AUROSTD_XSCALAR_TOLERANCE_INTEGER_);  //CO20191201

  template<class utype> bool _iszero(utype a,utype tol=(utype)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);  //CO20191201
  //bool iszero(bool x,bool tolerance); //CO20191201
  //bool iszero(char x,char tolerance); //CO20191201
  bool iszero(uint x,uint tolerance=(uint)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201
  bool iszero(int x,int tolerance=(int)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201
  bool iszero(long int x,long int tolerance=(long int)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201
  bool iszero(unsigned long int x,unsigned long int tolerance=(unsigned long int)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201
  bool iszero(long long int x,long long int tolerance=(long long int)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201
  bool iszero(unsigned long long int x,unsigned long long int tolerance=(unsigned long long int)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201
  bool iszero(float x,float tolerance=(float)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201
  bool iszero(double x,double tolerance=(double)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201
  bool iszero(long double x,long double tolerance=(long double)_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_); //CO20191201

  string dbl2frac(double a, bool sign_prefix=true); //DX20190724
  double frac2dbl(const string& str); //DX20200313
  void double2fraction(const double& input_double, int& numerator, int& denominator, double tol_diff=AUROSTD_IDENTITY_TOL, double tol_remainder=1e-2); //DX20210908
  int getNumeratorContinuedFractions(int& p, const int& n, vector<int>& fraction_sequence); //DX20210908
  int getDenominatorContinuedFractions(int& q, const int& n, vector<int>& fraction_sequence); //DX20210908

  template<class utype> utype fact(utype) __xprototype;
  template<class utype> utype factorial(utype) __xprototype;
  template<class utype> utype angle(utype,utype,utype,utype,utype,utype) __xprototype;
  template<class utype> utype angle(utype,utype,utype,utype) __xprototype;
  template<class utype> utype modulus(utype,utype,utype) __xprototype;
  template<class utype> utype modulus(utype,utype) __xprototype;
  // [OBSOLETE]  template<class utype> utype round(utype) __xprototype;
  // [OBSOLETE]  template<class utype> utype floor(utype) __xprototype;
  // [OBSOLETE]  template<class utype> utype ceil(utype) __xprototype;
  // [OBSOLETE]  template<class utype> utype trunc(utype) __xprototype;

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
}

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
#ifdef _AUROSTD_XCOMPLEX_
  bool _isfloat(xcomplex<float>) __xprototype;
  bool _isfloat(xcomplex<double>) __xprototype;
  bool _isfloat(xcomplex<long double>) __xprototype;
#endif
  bool isfloat(const string& in); //CO20180729
}

// ----------------------------------------------------------------------------
// _iscomplex  _iscomplex  _iscomplex  _iscomplex  _iscomplex
namespace aurostd {
  // namespace aurostd
  bool _iscomplex(bool) __xprototype;
  bool _iscomplex(char) __xprototype;
  bool _iscomplex(uint) __xprototype;
  bool _iscomplex(int) __xprototype;
  bool _iscomplex(long int) __xprototype;
  bool _iscomplex(unsigned long int) __xprototype;  //CO20191201
  bool _iscomplex(long long int) __xprototype;
  bool _iscomplex(unsigned long long int) __xprototype;
  bool _iscomplex(float) __xprototype;
  bool _iscomplex(double) __xprototype;
  bool _iscomplex(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  bool _iscomplex(xcomplex<float>) __xprototype;
  bool _iscomplex(xcomplex<double>) __xprototype;
  bool _iscomplex(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// _isreal  _isreal  _isreal  _isreal  _isreal
namespace aurostd {
  // namespace aurostd
  bool _isreal(bool) __xprototype;
  bool _isreal(char) __xprototype;
  bool _isreal(uint) __xprototype;
  bool _isreal(int) __xprototype;
  bool _isreal(long int) __xprototype;
  bool _isreal(unsigned long int) __xprototype; //CO20191201
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
}

// ----------------------------------------------------------------------------
// _isreal  _isreal  _isreal  _isreal  _isreal
namespace aurostd {
  // namespace aurostd
  bool _isreal(bool) __xprototype;
  bool _isreal(char) __xprototype;
  bool _isreal(uint) __xprototype;
  bool _isreal(int) __xprototype;
  bool _isreal(long int) __xprototype;
  bool _isreal(unsigned long int) __xprototype; //CO20191201
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
}

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
}

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
}

// ----------------------------------------------------------------------------
// _real _real _real _real _real _real _real _real
namespace aurostd {
  // namespace aurostd
  bool _real(bool) __xprototype;
  char _real(char) __xprototype;
  uint _real(uint) __xprototype;
  int _real(int) __xprototype;
  long int _real(long int) __xprototype;
  unsigned long int _real(unsigned long int) __xprototype;  //CO20191201
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
}

// ----------------------------------------------------------------------------
// sign sign sign sign sign sign sign sign
namespace aurostd {
  // namespace aurostd
  bool sign(bool) __xprototype;
  char sign(char) __xprototype;
  int sign(int) __xprototype;
  uint sign(uint) __xprototype;
  float sign(float) __xprototype;
  double sign(double) __xprototype;
  long int sign(long int) __xprototype;
  long long int sign(long long int) __xprototype;
  long double sign(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  xcomplex<float> sign(xcomplex<float>) __xprototype;
  xcomplex<double> sign(xcomplex<double>) __xprototype;
  xcomplex<long double> sign(xcomplex<long double>) __xprototype;
#endif
}

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
}

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
}

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
}

// ----------------------------------------------------------------------------
// sqrt  sqrt  sqrt  sqrt  sqrt
namespace aurostd {
  // namespace aurostd
  // SQRT(X)
  // int sqrt(int x) __xprototype;
  char sqrt(char x) __xprototype;
  int sqrt(int x) __xprototype;
  uint sqrt(uint x) __xprototype;
  float sqrt(float x) __xprototype;
  double sqrt(double x) __xprototype;
  long int sqrt(long int x) __xprototype;
  long long int sqrt(long long int x) __xprototype;
  unsigned long long int sqrt(unsigned long long int x) __xprototype;
  long double sqrt(long double x) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  //  float sqrt(xcomplex<float> x) __xprototype;
  //  double sqrt(xcomplex<double> x) __xprototype;
  long double sqrt(xcomplex<long double> x) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// exp  exp  exp  exp  exp
namespace aurostd {
  // namespace aurostd
  // EXP(X)
  char exp(char x) __xprototype; //SD20220512
  int exp(int x) __xprototype; //SD20220512
  uint exp(uint x) __xprototype; //SD20220512
  float exp(float x) __xprototype; //SD20220512
  double exp(double x) __xprototype; //SD20220512
  long int exp(long int x) __xprototype; //SD20220512
  long long int exp(long long int x) __xprototype; //SD20220512
  unsigned long long int exp(unsigned long long int x) __xprototype; //SD20220512
  long double exp(long double x) __xprototype; //SD20220512
}

// ----------------------------------------------------------------------------
// pow  pow  pow  pow  pow
namespace aurostd {
  // namespace aurostd
  // POW(X)
  char pow(char x,char d) __xprototype; //SD20220512
  int pow(int x,int d) __xprototype; //SD20220512
  uint pow(uint x, uint d) __xprototype; //SD20220512
  float pow(float x, float d) __xprototype; //SD20220512
  double pow(double x,double d) __xprototype; //SD20220512
  long int pow(long int x,long int d) __xprototype; //SD20220512
  long long int pow(long long int x,long long int d) __xprototype; //SD20220512
  unsigned long long int pow(unsigned long long int x,unsigned long long int d) __xprototype; //SD20220512
  long double pow(long double x,long double d) __xprototype; //SD20220512
}

//--------------------------------------------------------------- round
namespace aurostd {
  double round(double x, uint digits=0); //SD20220603
  int roundDouble(double doub, int multiple, bool up);  //CO20220624 (moved from chull)
  bool greaterEqualZero(double val);  //CO20220624 (moved from chull)
  bool lessEqualZero(double val); //CO20220624 (moved from chull)
  bool notPositive(double val,bool soft_cutoff,double tol);  //CO20220624 (moved from chull)
  bool notNegative(double val,bool soft_cutoff,double tol);  //CO20220624 (moved from chull)
  bool zeroWithinTol(double val,double tol); //CO20220624 (moved from chull)
  bool nonZeroWithinTol(double val,double tol);  //CO20220624 (moved from chull)
}

//--------------------------------------------------------------- isequal
namespace aurostd {
  // with const utype&
  template<class utype> bool identical(const utype&,const utype&,const utype&) __xprototype;
  template<class utype> bool identical(const utype&,const utype&) __xprototype;
  bool identical (const bool a,const bool b); //SD20220705
  bool identical (const char a,const char b); //SD20220705
  bool identical (const string& a,const string& b); //SD20220705
  template<class utype> bool isdifferent(const utype&,const utype&,const utype&) __xprototype;
  template<class utype> bool isdifferent(const utype&,const utype&) __xprototype;
  bool isdifferent (const bool a,const bool b); //SD20220705
  bool isdifferent (const char a,const char b); //SD20220705
  bool isdifferent (const string& a,const string& b); //SD20220705
  template<class utype> bool isequal(const utype&,const utype&,const utype&) __xprototype;
  template<class utype> bool isequal(const utype&,const utype&) __xprototype;
  bool isequal (const bool a,const bool b); //SD20220705
  bool isequal (const char a,const char b); //SD20220705
  bool isequal (const string& a,const string& b); //SD20220705
  // with utype
  //template<class utype> bool identical(utype,utype,utype) __xprototype;
  //template<class utype> bool identical(utype,utype) __xprototype;
  //template<class utype> bool isdifferent(utype,utype,utype) __xprototype;
  //template<class utype> bool isdifferent(utype,utype) __xprototype;
  //template<class utype> bool isequal(utype,utype,utype) __xprototype;
  //template<class utype> bool isequal(utype,utype) __xprototype;
}

//--------------------------------------------------------------- extra min/max // __XEXTRA_MINMAX_CPP
namespace aurostd {
  // namespace aurostd
  template<class utype> utype min(utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype min(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
  template<class utype> utype max(utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype,utype) __xprototype;  // __XEXTRA_MINMAX_CPP
}

namespace aurostd {
  template<class utype> utype _roundoff(const utype& x,utype tolerance);
  int roundoff(int x,int tolerance);
  long roundoff(long x,long tolerance);
  uint roundoff(uint x,uint tolerance);
  float roundoff(float x,float tolerance);
  double roundoff(double x,double tolerance);
  long long int roundoff(long long int x,long long int tolerance);
  unsigned long int roundoff(unsigned long int x,unsigned long int tolerance);
  unsigned long long int roundoff(unsigned long long int x,unsigned long long int tolerance);
  long double roundoff(long double x,long double tolerance);
}

namespace aurostd {
  int boundary_conditions_periodic(int lrows,int urows,int i);  //CO20190419 - taken from xvector BOUNDARY_CONDITIONS_PERIODIC
}

namespace aurostd {
  uint powint(uint x,uint exp); //CO20191201
  int powint(int x,uint exp); //CO20191201
  long int powint(long int x,uint exp); //CO20191201
  unsigned long int powint(unsigned long int x,uint exp); //CO20191201
  long long int powint(long long int x,uint exp); //CO20191201
  unsigned long long int powint(unsigned long long int x,uint exp); //CO20191201
}

//AS20200513 BEGIN
namespace aurostd{
  double FermiDirac(double E, double mu, double T);
}
//AS20200513 END

//CO20201111 BEGIN
namespace aurostd{
  template<class utype> utype nCk(utype n,utype k);
}
//CO20201111 END

//CO20201111 - BEGIN
namespace aurostd {
  bool isNaN(double d);
}
//CO20201111 - END

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XSCALAR_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

