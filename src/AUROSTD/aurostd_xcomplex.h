// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XCOMPLEX_H_
#define _AUROSTD_XCOMPLEX_H_

#define _AUROSTD_XCOMPLEX_

//*****************************************************************************
// -------------------------------------------------------------- class xcomplex
namespace aurostd {

  //[CO20200624 - OBSOLETE]class istream;
  //[CO20200624 - OBSOLETE]class ostream;

  template <class utype> class xcomplex;

  template <class utype>
    class xcomplex
    {
      public:
        xcomplex(const xcomplex& b);                          // constructor copy
        ~xcomplex();                                          // kill everything
        const xcomplex& operator=(const xcomplex &b);         // copy
        xcomplex(utype r=0,utype i=0): re(r),im(i) {;}        // constructors
        xcomplex& operator +=(const xcomplex&);
        xcomplex& operator -=(const xcomplex&);
        xcomplex& operator *=(const xcomplex&);
        xcomplex& operator /=(const xcomplex&);
        utype real() const { return re;}
        utype imag() const { return im;}
        //  private:
        utype re,im;
      private:
        //   friend std::ostream& operator<<(std::ostream &,const xcomplex&);    // print
        void free();
        // nothing because it is easier this way.
        //    utype re,im;
    };

  //ME20200107
  // ---------------------------------------------- identical
  template<class utype> bool identical(const xcomplex<utype>&, const xcomplex<utype>&, utype);
  template<class utype> bool isdifferent(const xcomplex<utype>&, const xcomplex<utype>&, utype);
  template<class utype> bool isdifferent(const xcomplex<utype>&, const xcomplex<utype>&);
  template<class utype> bool isequal(const xcomplex<utype>&, const xcomplex<utype>&, utype);
  template<class utype> bool isequal(const xcomplex<utype>&, const xcomplex<utype>&);

  // ---------------------------------------------- operator ostream << xcomplex
  //  template <class utype> ostream& operator <<(ostream&,const xcomplex<utype>&);
  template<class utype> std::ostream& operator <<(std::ostream& os,const xcomplex<utype>& x) __xprototype;
  // ---------------------------------------------- operator istream >> xcomplex
  template<class utype> std::istream& operator >>(std::istream&,xcomplex<utype>&);

  // ------------------------------------------------------------ function imag
  template <class utype> // removed inline
    utype imag(const xcomplex<utype>& x) __xprototype;

  // ------------------------------------------------------------ function real
  template <class utype> // removed inline
    utype real(const xcomplex<utype>& x) __xprototype;

  // ----------------------------------------------- operator xcomplex + xcomplex
  template <class utype> // removed inline
    xcomplex<utype> operator+(const xcomplex<utype>& x,const xcomplex<utype>& y) __xprototype;

  // ------------------------------------------------- operator xcomplex + utype
  template <class utype> // removed inline
    xcomplex<utype> operator+(const xcomplex<utype>& x,utype y) __xprototype;

  // ------------------------------------------------- operator utype + xcomplex
  template <class utype> // removed inline
    xcomplex<utype> operator+(utype x,const xcomplex<utype>& y) __xprototype;

  // ----------------------------------------------- operator xcomplex - xcomplex
  template <class utype> // removed inline
    xcomplex<utype> operator-(const xcomplex<utype>& x,const xcomplex<utype>& y) __xprototype;
  // ------------------------------------------------- operator xcomplex - utype
  template <class utype> // removed inline
    xcomplex<utype> operator-(const xcomplex<utype>& x,utype y) __xprototype;

  // ------------------------------------------------- operator utype - xcomplex
  template <class utype> // removed inline
    xcomplex<utype> operator-(utype x,const xcomplex<utype>& y) __xprototype;
  // ------------------------------------------------- operator xcomplex * utype
  template <class utype> // removed inline
    xcomplex<utype> operator*(const xcomplex<utype>& x,const xcomplex<utype>& y) __xprototype;

  // ------------------------------------------------- operator xcomplex * utype
  template <class utype> // removed inline
    xcomplex<utype> operator*(const xcomplex<utype>& x,utype y) __xprototype;

  // ------------------------------------------------- operator utype * xcomplex
  template <class utype> // removed inline
    xcomplex<utype> operator*(utype x,const xcomplex<utype>& y) __xprototype;

  // ----------------------------------------------- operator xcomplex / xcomplex
  template <class utype>
    xcomplex<utype> operator /(const xcomplex<utype>&,const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------- operator utype / xcomplex
  template <class utype> xcomplex<utype>
    operator /(utype,const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------- operator xcomplex / utype
  template <class utype> xcomplex<utype>
    operator /(const xcomplex<utype>& x,utype y) __xprototype;

  // -------------------------------------------------------- operator +xcomplex
  template <class utype> // removed inline
    xcomplex<utype> operator+(const xcomplex<utype>& x) __xprototype;

  // -------------------------------------------------------- operator -xcomplex
  template <class utype> // removed inline
    xcomplex<utype> operator-(const xcomplex<utype>& x) __xprototype;

  // ---------------------------------------------- operator xcomplex == xcomplex
  template <class utype> // removed inline
    bool operator ==(const xcomplex<utype>& x,const xcomplex<utype>& y) __xprototype;

  // ------------------------------------------------ operator xcomplex == utype
  template <class utype> // removed inline
    bool operator ==(const xcomplex<utype>& x,utype y)  __xprototype;

  // ------------------------------------------------ operator utype == xcomplex
  template <class utype> // removed inline
    bool operator ==(utype x,const xcomplex<utype>& y)  __xprototype;

  // ---------------------------------------------- operator xcomplex != xcomplex
  template <class utype> // removed inline
    bool operator !=(const xcomplex<utype>& x,const xcomplex<utype>& y)  __xprototype;

  // ------------------------------------------------ operator xcomplex != utype
  template <class utype> // removed inline
    bool  operator !=(const xcomplex<utype>& x,utype y)  __xprototype;

  // ------------------------------------------------ operator utype != xcomplex
  template <class utype> // removed inline
    bool  operator !=(utype x,const xcomplex<utype>& y)  __xprototype;

  // ------------------------------------------------------------- function abs
  template <class utype> // removed inline
    utype abs(const xcomplex<utype>& x) __xprototype;

  // ------------------------------------------------------------- function arg
  template <class utype> // removed inline
    utype  arg(const xcomplex<utype>& x) __xprototype;

  // ----------------------------------------------------------- function polar
  template <class utype> // removed inline
    xcomplex<utype>  polar(utype r,utype t) __xprototype;

  // ------------------------------------------------------------ function conj
  template <class utype> // removed inline
    xcomplex<utype>  conj(const xcomplex<utype>& x) __xprototype;
  template <class utype> //ME20180904
    utype  conj(const utype& x) __xprototype;
  // ------------------------------------------------------------ function norm
  template <class utype> // removed inline
    utype norm(const xcomplex<utype>& x) __xprototype;
  // ------------------------------------------------------------- function cos
  template <class utype> xcomplex<utype>
    cos(const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------------------ function cosh
  template <class utype> xcomplex<utype>
    cosh(const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------------------- function exp
  template <class utype> xcomplex<utype>
    exp(const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------------------- function log
  template <class utype> xcomplex<utype>
    log(const xcomplex<utype>&) __xprototype;

  // --------------------------------------------- function xcomplex pow xcomplex
  template <class utype> xcomplex<utype>
    pow(const xcomplex<utype>&,const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------- function xcomplex pow int
  template <class utype> xcomplex<utype>
    pow(const xcomplex<utype>&,int) __xprototype;

  // ----------------------------------------------- function xcomplex pow utype
  template <class utype> xcomplex<utype>
    pow(const xcomplex<utype>&,utype) __xprototype;

  // ----------------------------------------------- function utype pow xcomplex
  template <class utype> xcomplex<utype>
    pow(utype,const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------------------- function sin
  template <class utype> xcomplex<utype>
    sin(const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------------------ function sinh
  template <class utype> xcomplex<utype>
    sinh(const xcomplex<utype>&) __xprototype;

  // ------------------------------------------------------------ function sqrt
  template <class utype> xcomplex<utype>
    sqrt(const xcomplex<utype>&) __xprototype;

  //ME20180907
  // ------------------------------------------------------------ function magnitude square
  template <class utype> utype
    magsqr(const xcomplex<utype>&) __xprototype;

}

// // ---------------------------------------------------------------------------- 
// // ------------------------------------------ specialization for xcomplex<float> 

// namespace aurostd { 
//   // template <class utype>  
//   template <> class xcomplex<float> 
//   { 
//   public: 
//     xcomplex(float r=0,float i=0): re(r),im(i) { } 
//     explicit xcomplex (const xcomplex<double>& r); 
//     explicit xcomplex (const xcomplex<long double>& r); 

//     xcomplex& operator+=(const xcomplex& r); 
//     xcomplex& operator-=(const xcomplex& r); 
//     xcomplex& operator*=(const xcomplex& r); 
//     xcomplex& operator/=(const xcomplex& r); 

//     float real() const { return re;} 
//     float imag() const { return im;} 
//   private: 
//     float re,im; 

// #ifndef __STRICT_ANSI__ 
//     friend // removed inline  
//     xcomplex operator+(const xcomplex& x,float y) 
//     { return operator+<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator+(float x,const xcomplex& y) 
//     { return operator+<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator-(const xcomplex& x,float y) 
//     { return operator-<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator-(float x,const xcomplex& y) 
//     { return operator-<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator*(const xcomplex& x,float y) 
//     { return operator*<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator*(float x,const xcomplex& y) 
//     { return operator*<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator /(const xcomplex& x,float y) 
//     { return operator/<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator /(float x,const xcomplex& y) 
//     { return operator/<>(x,y);} 
//     friend // removed inline  
//     bool operator ==(const xcomplex& x,float y) 
//     { return operator==<>(x,y);} 
//     friend // removed inline  
//     bool operator ==(float x,const xcomplex& y) 
//     { return operator==<>(x,y);} 
//     friend // removed inline  
//     bool operator !=(const xcomplex& x,float y) 
//     { return operator!=<>(x,y);} 
//     friend // removed inline  
//     bool operator !=(float x,const xcomplex& y) 
//     { return operator!=<>(x,y);} 
// #endif /\* __STRICT_ANSI__ *\/ 
//   }; 
// } 

// // ---------------------------------------------------------------------------- 
// // ----------------------------------------- specialization for xcomplex<double> 

// namespace aurostd { 
//   template <> class xcomplex<double> 
//   { 
//   public: 
//     xcomplex(double r=0,double i=0): re(r),im(i) { } 
//     xcomplex(const xcomplex<float>& r): re(r.real()),im(r.imag()) { } 
//     explicit xcomplex (const xcomplex<long double>& r); 

//     xcomplex& operator+=(const xcomplex& r); 
//     xcomplex& operator-=(const xcomplex& r); 
//     xcomplex& operator*=(const xcomplex& r); 
//     xcomplex& operator/=(const xcomplex& r); 

//     double real() const { return re;} 
//     double imag() const { return im;} 
//   private: 
//     double re,im; 

// #ifndef __STRICT_ANSI__ 
//     friend // removed inline  
//     xcomplex operator+(const xcomplex& x,double y) 
//     { return operator+<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator+(double x,const xcomplex& y) 
//     { return operator+<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator-(const xcomplex& x,double y) 
//     { return operator-<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator-(double x,const xcomplex& y) 
//     { return operator-<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator*(const xcomplex& x,double y) 
//     { return operator*<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator*(double x,const xcomplex& y) 
//     { return operator*<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator /(const xcomplex& x,double y) 
//     { return operator/<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator /(double x,const xcomplex& y) 
//     { return operator/<>(x,y);} 
//     friend // removed inline  
//     bool operator ==(const xcomplex& x,double y) 
//     { return operator==<>(x,y);} 
//     friend // removed inline  
//     bool operator ==(double x,const xcomplex& y) 
//     { return operator==<>(x,y);} 
//     friend // removed inline  
//     bool operator !=(const xcomplex& x,double y) 
//     { return operator!=<>(x,y);} 
//     friend // removed inline  
//     bool operator !=(double x,const xcomplex& y) 
//     { return operator!=<>(x,y);} 
// #endif /\* __STRICT_ANSI__ *\/ 
//   }; 

//   // removed inline  
//   xcomplex<float>::xcomplex(const xcomplex<double>& r) 
//     : re(r.real()),im(r.imag()) 
//   { } 
// } 

// // ---------------------------------------------------------------------------- 
// // ------------------------------------ specialization for xcomplex<long double> 

// namespace aurostd { 
//   template <> class xcomplex<long double> 
//   { 
//   public: 
//     xcomplex(long double r=0,long double i=0): re(r),im(i) { } 
//     xcomplex(const xcomplex<float>& r): re(r.real()),im(r.imag()) { } 
//     xcomplex(const xcomplex<double>& r): re(r.real()),im(r.imag()) { } 

//     xcomplex& operator+=(const xcomplex& r); 
//     xcomplex& operator-=(const xcomplex& r); 
//     xcomplex& operator*=(const xcomplex& r); 
//     xcomplex& operator/=(const xcomplex& r); 

//     long double real() const { return re;} 
//     long double imag() const { return im;} 
//   private: 
//     long double re,im; 

// #ifndef __STRICT_ANSI__ 
//     friend // removed inline  
//     xcomplex operator+(const xcomplex& x,long double y) 
//     { return operator+<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator+(long double x,const xcomplex& y) 
//     { return operator+<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator-(const xcomplex& x,long double y) 
//     { return operator-<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator-(long double x,const xcomplex& y) 
//     { return operator-<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator*(const xcomplex& x,long double y) 
//     { return operator*<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator*(long double x,const xcomplex& y) 
//     { return operator*<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator /(const xcomplex& x,long double y) 
//     { return operator/<>(x,y);} 
//     friend // removed inline  
//     xcomplex operator /(long double x,const xcomplex& y) 
//     { return operator/<>(x,y);} 
//     friend // removed inline  
//     bool operator ==(const xcomplex& x,long double y) 
//     { return operator==<>(x,y);} 
//     friend // removed inline  
//     bool operator ==(long double x,const xcomplex& y) 
//     { return operator==<>(x,y);} 
//     friend // removed inline  
//     bool operator !=(const xcomplex& x,long double y) 
//     { return operator!=<>(x,y);} 
//     friend // removed inline  
//     bool operator !=(long double x,const xcomplex& y) 
//     { return operator!=<>(x,y);} 
// #endif /\* __STRICT_ANSI__ *\/ 
//   }; 

//   // removed inline  
//   xcomplex<float>::xcomplex(const xcomplex<long double>& r) 
//     : re(r.real()),im(r.imag()) 
//   { } 

//   // removed inline  
//   xcomplex<double>::xcomplex(const xcomplex<long double>& r) 
//     : re(r.real()),im(r.imag()) 
//   { } 
// } // extern "C++" 

// // ---------------------------------------------------------------------------- 
// // ---------------------------------------------------------------------------- 
// // ---------------------------------------------------------------------------- 

// //typedef xcomplex<float> float_xcomplex; 
// //typedef xcomplex<double> double_xcomplex; 
// //typedef xcomplex<long double> long_double_xcomplex; 

#endif // _AUROSTD_XCOMPLEX_H_


// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
