// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2012

#ifndef _AUROSTD_XRANDOM_CPP_
#define _AUROSTD_XRANDOM_CPP_

#ifndef XXEND
#define XXEND 1
#endif

#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

// ----------------------------------------------------------------------------
// --------------------------------------------------------------- implementation
#define RAN3
static long int _idum[1];                                    // to start ..

namespace aurostd {
  // namespace aurostd
  long int _random_initialize(long int _num) {
    (*_idum)=-std::labs(_num)-1969;
    // (*_idum)=-abs(_num)-1969;
    return (*_idum);
  }

  long int _random_initialize(void) {
    //    (*_idum)=1969;
    //    (*_idum)=1969-aurostd::execute2int("date | sum")+aurostd::execute2int("ls /tmp | sum");

    timeval tim;
    gettimeofday(&tim, NULL);
    //   double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
    //   cerr << tim.tv_usec << endl;
    // getpid
    int PID=getpid();
    int TID=aurostd::getTID(); //CO20200502 - threadID
    long int seed=tim.tv_usec+tim.tv_sec+PID+TID; //CO20200502 - threadID
    srand(seed);      // for std:: library things
    (*_idum)=seed;    // for aurostd:: library things
    //    cerr << PID << endl;
    //    cerr << TID << endl;  //CO20200502 - threadID
    return (*_idum);
  }

  // --------------------------------------------------------- uniform deviations
  // --------------------------------------------------------------------- ran0()

#define _AUROSTD_XRANDOM_RAN0_IA 16807
#define _AUROSTD_XRANDOM_RAN0_IM 2147483647
#define _AUROSTD_XRANDOM_RAN0_AM (1.0/_AUROSTD_XRANDOM_RAN0_IM)
#define _AUROSTD_XRANDOM_RAN0_IQ 127773
#define _AUROSTD_XRANDOM_RAN0_IR 2836
#define _AUROSTD_XRANDOM_RAN0_MASK 123459876

  double ran0(void) { // 108.3 ns
    long int k;
    float ans;

    *_idum ^= _AUROSTD_XRANDOM_RAN0_MASK;
    k=(*_idum)/_AUROSTD_XRANDOM_RAN0_IQ;
    *_idum=_AUROSTD_XRANDOM_RAN0_IA*(*_idum-k*_AUROSTD_XRANDOM_RAN0_IQ)-_AUROSTD_XRANDOM_RAN0_IR*k;
    if(*_idum < 0) *_idum += _AUROSTD_XRANDOM_RAN0_IM;
    ans=_AUROSTD_XRANDOM_RAN0_AM*(*_idum);
    *_idum ^= _AUROSTD_XRANDOM_RAN0_MASK;
    return ans;
  }

#undef _AUROSTD_XRANDOM_RAN0_IA
#undef _AUROSTD_XRANDOM_RAN0_IM
#undef _AUROSTD_XRANDOM_RAN0_AM
#undef _AUROSTD_XRANDOM_RAN0_IQ
#undef _AUROSTD_XRANDOM_RAN0_IR
#undef _AUROSTD_XRANDOM_RAN0_MASK

  // --------------------------------------------------------------------- ran1()
#define _AUROSTD_XRANDOM_RAN1_IA        16807
#define _AUROSTD_XRANDOM_RAN1_IM        2147483647
#define _AUROSTD_XRANDOM_RAN1_AM        (1.0/_AUROSTD_XRANDOM_RAN1_IM)
#define _AUROSTD_XRANDOM_RAN1_IQ        127773
#define _AUROSTD_XRANDOM_RAN1_IR        2836
#define _AUROSTD_XRANDOM_RAN1_NTAB      32
#define _AUROSTD_XRANDOM_RAN1_EPS       1.2e-7
#define _AUROSTD_XRANDOM_RAN1_RNMX      (1.0-_AUROSTD_XRANDOM_RAN1_EPS)

  double ran1(void) { // 118.80 ns
    // "minimal" random number generator of Park and Miller with Bays-Durham
    // shuffle and added safeguards. Return a uniform random deviate between
    // 0.0 and 1.0 (exclusive of the endpoint values:
    //  _AUROSTD_XRANDOM_RAN1_EPS and _AUROSTD_XRANDOM_RAN1_RNMX)

    int j;
    long int k;
    static long int iy=0;
    static long int iv[_AUROSTD_XRANDOM_RAN1_NTAB];
    double temp;
    if(*_idum<=0 || !iy) {
      if(-(*_idum)<1)
        *_idum=1;
      else
        *_idum=-(*_idum);
      for(j=_AUROSTD_XRANDOM_RAN1_NTAB+7;j>=0;j--) {
        k=(*_idum)/_AUROSTD_XRANDOM_RAN1_IQ;
        *_idum=_AUROSTD_XRANDOM_RAN1_IA*(*_idum-k*_AUROSTD_XRANDOM_RAN1_IQ)-_AUROSTD_XRANDOM_RAN1_IR*k;
        if(*_idum<0)
          *_idum+=_AUROSTD_XRANDOM_RAN1_IM;
        if(j<_AUROSTD_XRANDOM_RAN1_NTAB)
          iv[j]=*_idum;
      }
      iy=iv[0];
    }
    k=(*_idum)/_AUROSTD_XRANDOM_RAN1_IQ;
    *_idum=_AUROSTD_XRANDOM_RAN1_IA*(*_idum-k*_AUROSTD_XRANDOM_RAN1_IQ)-_AUROSTD_XRANDOM_RAN1_IR*k;
    if(*_idum<0)
      *_idum+=_AUROSTD_XRANDOM_RAN1_IM;
    j=iy/(1+(_AUROSTD_XRANDOM_RAN1_IM-1)/_AUROSTD_XRANDOM_RAN1_NTAB);
    iy=iv[j];
    iv[j] = *_idum;
    if((temp=_AUROSTD_XRANDOM_RAN1_AM*iy)>_AUROSTD_XRANDOM_RAN1_RNMX)
      return _AUROSTD_XRANDOM_RAN1_RNMX;
    else
      return temp;
  }

#undef _AUROSTD_XRANDOM_RAN1_IA
#undef _AUROSTD_XRANDOM_RAN1_IM
#undef _AUROSTD_XRANDOM_RAN1_AM
#undef _AUROSTD_XRANDOM_RAN1_IQ
#undef _AUROSTD_XRANDOM_RAN1_IR
#undef _AUROSTD_XRANDOM_RAN1_NTAB
#undef _AUROSTD_XRANDOM_RAN1_EPS
#undef _AUROSTD_XRANDOM_RAN1_RNMX

  // --------------------------------------------------------------------- ran2()

#define _AUROSTD_XRANDOM_RAN2_IM1 2147483563
#define _AUROSTD_XRANDOM_RAN2_IM2 2147483399
#define _AUROSTD_XRANDOM_RAN2_AM (1.0/_AUROSTD_XRANDOM_RAN2_IM1)
#define _AUROSTD_XRANDOM_RAN2_IMM1 (_AUROSTD_XRANDOM_RAN2_IM1-1)
#define _AUROSTD_XRANDOM_RAN2_IA1 40014
#define _AUROSTD_XRANDOM_RAN2_IA2 40692
#define _AUROSTD_XRANDOM_RAN2_IQ1 53668
#define _AUROSTD_XRANDOM_RAN2_IQ2 52774
#define _AUROSTD_XRANDOM_RAN2_IR1 12211
#define _AUROSTD_XRANDOM_RAN2_IR2 3791
#define _AUROSTD_XRANDOM_RAN2_NTAB 32
#define _AUROSTD_XRANDOM_RAN2_NDIV (1+_AUROSTD_XRANDOM_RAN2_IMM1/_AUROSTD_XRANDOM_RAN2_NTAB)
#define _AUROSTD_XRANDOM_RAN2_EPS 1.2e-7
#define _AUROSTD_XRANDOM_RAN2_RNMX (1.0-_AUROSTD_XRANDOM_RAN2_EPS)

  double ran2(void) { // 215.0 ns
    int j;
    long int k;
    static long int _idum2=123456789;
    static long int iy=0;
    static long int iv[_AUROSTD_XRANDOM_RAN2_NTAB];
    float temp;

    if(*_idum <= 0) {
      if(-(*_idum) < 1) *_idum=1;
      else *_idum = -(*_idum);
      _idum2=(*_idum);
      for (j=_AUROSTD_XRANDOM_RAN2_NTAB+7;j>=0;j--) {
        k=(*_idum)/_AUROSTD_XRANDOM_RAN2_IQ1;
        *_idum=_AUROSTD_XRANDOM_RAN2_IA1*(*_idum-k*_AUROSTD_XRANDOM_RAN2_IQ1)-k*_AUROSTD_XRANDOM_RAN2_IR1;
        if(*_idum < 0) *_idum += _AUROSTD_XRANDOM_RAN2_IM1;
        if(j < _AUROSTD_XRANDOM_RAN2_NTAB) iv[j] = *_idum;
      }
      iy=iv[0];
    }
    k=(*_idum)/_AUROSTD_XRANDOM_RAN2_IQ1;
    *_idum=_AUROSTD_XRANDOM_RAN2_IA1*(*_idum-k*_AUROSTD_XRANDOM_RAN2_IQ1)-k*_AUROSTD_XRANDOM_RAN2_IR1;
    if(*_idum < 0) *_idum += _AUROSTD_XRANDOM_RAN2_IM1;
    k=_idum2/_AUROSTD_XRANDOM_RAN2_IQ2;
    _idum2=_AUROSTD_XRANDOM_RAN2_IA2*(_idum2-k*_AUROSTD_XRANDOM_RAN2_IQ2)-k*_AUROSTD_XRANDOM_RAN2_IR2;
    if(_idum2 < 0) _idum2 += _AUROSTD_XRANDOM_RAN2_IM2;
    j=iy/_AUROSTD_XRANDOM_RAN2_NDIV;
    iy=iv[j]-_idum2;
    iv[j] = *_idum;
    if(iy < 1) iy += _AUROSTD_XRANDOM_RAN2_IMM1;
    if((temp=_AUROSTD_XRANDOM_RAN2_AM*iy) > _AUROSTD_XRANDOM_RAN2_RNMX)
      return _AUROSTD_XRANDOM_RAN2_RNMX;
    else return temp;
  }

#undef _AUROSTD_XRANDOM_RAN2_IM1
#undef _AUROSTD_XRANDOM_RAN2_IM2
#undef _AUROSTD_XRANDOM_RAN2_AM
#undef _AUROSTD_XRANDOM_RAN2_IMM1
#undef _AUROSTD_XRANDOM_RAN2_IA1
#undef _AUROSTD_XRANDOM_RAN2_IA2
#undef _AUROSTD_XRANDOM_RAN2_IQ1
#undef _AUROSTD_XRANDOM_RAN2_IQ2
#undef _AUROSTD_XRANDOM_RAN2_IR1
#undef _AUROSTD_XRANDOM_RAN2_IR2
#undef _AUROSTD_XRANDOM_RAN2_NTAB
#undef _AUROSTD_XRANDOM_RAN2_NDIV
#undef _AUROSTD_XRANDOM_RAN2_EPS
#undef _AUROSTD_XRANDOM_RAN2_RNMX

  // --------------------------------------------------------------------- ran3()

#define _AUROSTD_XRANDOM_RAN3_MBIG 1000000000
#define _AUROSTD_XRANDOM_RAN3_MSEED 161803398
#define _AUROSTD_XRANDOM_RAN3_MZ 0
#define _AUROSTD_XRANDOM_RAN3_FAC (1.0/_AUROSTD_XRANDOM_RAN3_MBIG)

  double ran3(void) { // 83.60 ns
    static int inext,inextp;
    static long int ma[56];
    static int iff=0;
    long int mj,mk;
    int i,ii,k;

    if(*_idum < 0 || iff == 0) {
      iff=1;
      mj=_AUROSTD_XRANDOM_RAN3_MSEED-(*_idum < 0 ? -*_idum : *_idum);
      mj %= _AUROSTD_XRANDOM_RAN3_MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
        ii=(21*i) % 55;
        ma[ii]=mk;
        mk=mj-mk;
        if(mk < _AUROSTD_XRANDOM_RAN3_MZ) mk += _AUROSTD_XRANDOM_RAN3_MBIG;
        mj=ma[ii];
      }
      for (k=1;k<=4;k++)
        for (i=1;i<=55;i++) {
          ma[i] -= ma[1+(i+30) % 55];
          if(ma[i] < _AUROSTD_XRANDOM_RAN3_MZ) ma[i] += _AUROSTD_XRANDOM_RAN3_MBIG;
        }
      inext=0;
      inextp=31;
      *_idum=1;
    }
    if(++inext == 56) inext=1;
    if(++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if(mj < _AUROSTD_XRANDOM_RAN3_MZ) mj += _AUROSTD_XRANDOM_RAN3_MBIG;
    ma[inext]=mj;
    return mj*_AUROSTD_XRANDOM_RAN3_FAC;
  }
#undef _AUROSTD_XRANDOM_RAN3_MBIG
#undef _AUROSTD_XRANDOM_RAN3_MSEED
#undef _AUROSTD_XRANDOM_RAN3_MZ
#undef _AUROSTD_XRANDOM_RAN3_FAC

  // -------------------------------------------------------------- generic ran()

#ifdef RAN0
  inline double ran(void) { return ran0(); }
#else
#ifdef RAN1
  inline double ran(void) { return ran1(); }
#else
#ifdef RAN2
  inline double ran(void) { return ran2(); }
#else
#ifdef RAN3
  inline double ran(void) { return ran3(); }
#else
  inline double ran(void) { return ran0(); }
#endif
#endif
#endif
#endif

  template<class utype> inline utype
    uniform(const utype& x) {
      // returns a value from 0 to x, using ran() algorithm  ....
      return (utype) x*((utype) ran());
    }

  template<class utype> inline utype
    uniform(const utype& x,const utype& y) {
      // returns a value from x to y, using ran() algorithm  ....
      return (utype) x+(y-x)*((utype) ran());
    }


  //   inline int uniform(const int& x) { return (int) x*((int) ran());}
  //   inline int uniform(const int& x,const int& y) { return (int) x+(y-x)*((int) ran());}
  //   inline float uniform(const float& x) { return (float) x*((float) ran());}
  //   inline float uniform(const float& x,const float& y) { return (float) x+(y-x)*((float) ran());}
  //   inline double uniform(const double& x) { return (double) x*((double) ran());}
  //   inline double uniform(const double& x,const double& y) { return (double) x+(y-x)*((double) ran());}

  // --------------------------------------------------------- gaussian deviation

  double gaussian(void) {                                // Normal mean=0 sigma=1
    // returns a normally distribuited deviate with zero mean and unit
    // variance, using ran() as source of uniform deviates
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

    if(iset==0) {
      do {                                               // we need two
        v1=2.0*ran()-1.0;                               // uniforms
        v2=2.0*ran()-1.0;                               // variabiles
        rsq=v1*v1+v2*v2;                                 // in the unit
      } while (rsq>=1.0 || rsq == 0.0);                  // circle ..
      fac=std::sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;                                            //  set flag
      return v2*fac;
    } else {
      iset=0;                                            // clear flag
      return gset;
    }
  }

  template<class utype> inline utype                    // gaussian(mean=0,sigma)
    gaussian(const utype& sigma) {
      // returns a gaussian deviate with mean=0, sigma. using gaussian()
      //  return (utype) sqrt(sigma)*gaussian();
      return (utype) sigma*((utype) gaussian());
    }

  template<class utype> inline utype                         // gaussian(m,sigma)
    gaussian(const utype& mean,const utype& sigma) {
      // returns a gaussian deviate with mean=mean, sigma. using gaussian()
      //  return (utype) mean+sqrt(sigma)*gaussian();
      return (utype) mean+sigma*((utype) gaussian());
    }

  // ------------------------------------------------------ exponential deviation

  inline double expdev(void) {                            // exponential lambda=1
    // returns a exponentially distribuited deviate with lambda=1
    double u;  
    do u=ran();
    while (u == 0.0);
    return -log(u);
  }

  template<class utype> inline utype                        // exponential lambda
    expdev(const utype& lambda) {
      // returns a exponentially distribuited deviate with lambda
      double u;  
      do u=ran();
      while (u == 0.0);
      return (utype) -lambda*((utype) log(u));
    }

  // ---------------------------------------------------------- laplace deviation

  inline double laplacedev(void) {                          // laplace mean=0 a=1
    // returns a laplace distribuited deviate p(x)=1/2a*exp(-|x-lambda|/a)
    // sigma^2=2a*a mean=lambda
    double u;
    do u=ran();
    while (u == 0.0);
    return u<0.5 ? log(2.0*u):-log(2.0*(1.0-u));
  }

  template<class utype> inline utype                          // laplace mean=0 a
    laplacedev(const utype& a) {
      // returns a laplace distribuited deviate p(x)=1/2a*exp(-|x-lambda|/a)
      // sigma^2=2a*a mean=lambda
      double u;
      do u=ran();
      while (u == 0.0);
      return u<0.5 ? a*((utype) log(2.0*u)):-a*((utype) log(2.0*(1.0-u)));
    }

  template<class utype> inline utype                            // laplace mean a
    laplacedev(const utype& a,const utype& lambda) {
      // returns a laplace distribuited deviate p(x)=1/2a*exp(-|x-lambda|/a)
      // sigma^2=2a*a mean=lambda
      double u;
      do u=ran();
      while (u == 0.0);
      return u<0.5 ? a*((utype) log(2.0*u))+lambda:-a*((utype) log(2.0*(1.0-u)))+lambda;
    }

  // ---------------------------------------------------------- poisson deviation
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#ifdef _AUROSTD_XVECTOR_CPP_

namespace aurostd {
  // namespace aurostd
  template<class utype> xvector<utype>              // uniform xvector<utype> [x,y]
    uniform(const xvector<utype>& a,const utype& x,const utype& y)  {
      for(int i=a.lrows;i<=a.urows;i++)
        a[i]=(utype) x+(y-x)*uniform(utype(1.0));
      return a;
    }
 template xvector<double> uniform(const xvector<double>& a,const double& x,const double& y);

  template<class utype> xvector<utype>              // uniform xvector<utype> [0,y]
    uniform(const xvector<utype>& a,const utype& y)  {
      return uniform(a,utype(0.0),utype(y));
    }

  template<class utype> xvector<utype>              // uniform xvector<utype> [0,1]
    uniform(const xvector<utype>& a)  {
      return uniform(a,utype(0.0),utype(1.0));
    }


  template<class utype> xvector<xcomplex<utype> >     // xcomplex xvector of uniforms
    uniform(const xvector<xcomplex<utype> >& a,     // real,imag in ([x1,y1],[x2,y2])
        const utype & x1,const utype & y1,
        const utype & x2,const utype & y2) {
      for(int i=a.lrows;i<=a.urows;i++)
        a[i]=xcomplex<utype>(x1+(y1-x1)*uniform(utype(1.0)),
            x2+(y2-x2)*uniform(utype(1.0)));
      return a;
    }

  template<class utype> xvector<xcomplex<utype> >     // xcomplex xvector of uniforms
    uniform(const xvector<xcomplex<utype> >& a,       // real,imag in ([0,x],[0,y])
        const utype & x,const utype & y) {
      return uniform(a,utype(0.0),x,utype(0.0),y);
    }

  template<class utype> xvector<xcomplex<utype> >     // xcomplex xvector of uniforms
    uniform(const xvector<xcomplex<utype> >& a,       // real,imag in ([0,y],[0,y])
        const utype & y) {
      return uniform(a,utype(0.0),y,utype(0.0),y);
    }

  template<class utype> xvector<xcomplex<utype> >     // xcomplex xvector of uniforms
    uniform(const xvector<xcomplex<utype> >& a) {     // real,imag in ([0,1],[0,1])
      return uniform(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

// ----------------------------------------------------- gaussian distributions

namespace aurostd {
  // namespace aurostd
  template<class utype> xvector<utype>         // gaussian xvector<utype> (m,sigma)
    gaussian(const xvector<utype>& a,const utype& m,const utype& sigma)  {
      for(int i=a.lrows;i<=a.urows;i++)
        //      a[i]=(utype) m+sqrt(sigma)*gaussian();
        a[i]=(utype) m+sigma*gaussian();
      return a;
    }

  template<class utype> xvector<utype>         // gaussian xvector<utype> (m,sigma)
    gaussian(const xvector<utype>& a,
        const xvector<utype>& m,
        const xvector<utype>& sigma)  {
      for(int i=a.lrows;i<=a.urows;i++)
        //      a[i]=(utype) m[i]+sqrt(sigma[i])*gaussian();
        a[i]=(utype) m[i]+sigma[i]*gaussian();
      return a;
    }

  template<class utype> xvector<utype>         // gaussian xvector<utype> (0,sigma)
    gaussian(const xvector<utype>& a,const utype& sigma)  {
      return gaussian(a,utype(0.0),utype(sigma));
    }

  template<class utype> xvector<utype>             // gaussian xvector<utype> (0,1)
    gaussian(const xvector<utype>& a)  {
      return gaussian(a,utype(0.0),utype(1.0));
    }


  template<class utype> xvector<xcomplex<utype> >    // xcomplex xvector of gaussians
    gaussian(const xvector<xcomplex<utype> >& a, // real,imag (m1,sigma1),(m2,sigma2)
        const utype & m1,const utype & sigma1,
        const utype & m2,const utype & sigma2) {
      for(int i=a.lrows;i<=a.urows;i++)
        // a[i]=xcomplex<utype>(m1+sqrt(sigma1)*gaussian(),
        //                       m2+sqrt(sigma2)*gaussian());
        a[i]=xcomplex<utype>(m1+sigma1*gaussian(),
            m2+sigma2*gaussian());
      return a;
    }

  template<class utype> xvector<xcomplex<utype> >    // xcomplex xvector of gaussians
    gaussian(const xvector<xcomplex<utype> >& a,    // real,imag in (0,m),(0,sigma)
        const utype & m,const utype & sigma) {
      return gaussian(a,utype(0.0),m,utype(0.0),sigma);
    }

  template<class utype> xvector<xcomplex<utype> >    // xcomplex xvector of gaussians
    gaussian(const xvector<xcomplex<utype> >& a,// real,imag in (0,sigma),(0,sigma)
        const utype & sigma) {
      return gaussian(a,utype(0.0),sigma,utype(0.0),sigma);
    }

  template<class utype> xvector<xcomplex<utype> >    // xcomplex xvector of gaussians
    gaussian(const xvector<xcomplex<utype> >& a) {      // real,imag in (0,1),(0,1)
      return gaussian(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

#endif

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
#ifdef _AUROSTD_XMATRIX_CPP_

namespace aurostd {
  // namespace aurostd
  template<class utype> xmatrix<utype>              // uniform xmatrix<utype> [x,y]
    uniform(const xmatrix<utype>& a,const utype& x,const utype& y)  {
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          a[i][j]=(utype) x+(y-x)*uniform(utype(1.0));
      return a;
    }

  template<class utype> xmatrix<utype>              // uniform xmatrix<utype> [0,y]
    uniform(const xmatrix<utype>& a,const utype& y)  {
      return uniform(a,utype(0.0),utype(y));
    }

  template<class utype> xmatrix<utype>              // uniform xmatrix<utype> [0,1]
    uniform(const xmatrix<utype>& a)  {
      return uniform(a,utype(0.0),utype(1.0));
    }


  template<class utype> xmatrix<xcomplex<utype> >     // xcomplex xmatrix of uniforms
    uniform(const xmatrix<xcomplex<utype> >& a,   // real,imag in ([x1,y1],[x2,y2])
        const utype & x1,const utype & y1,
        const utype & x2,const utype & y2) {
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          a[i][j]=xcomplex<utype>(x1+(y1-x1)*uniform(utype(1.0)),
              x2+(y2-x2)*uniform(utype(1.0)));
      return a;
    }

  template<class utype> xmatrix<xcomplex<utype> >     // xcomplex xmatrix of uniforms
    uniform(const xmatrix<xcomplex<utype> >& a,       // real,imag in ([0,x],[0,y])
        const utype & x,const utype & y) {
      return uniform(a,utype(0.0),x,utype(0.0),y);
    }

  template<class utype> xmatrix<xcomplex<utype> >     // xcomplex xmatrix of uniforms
    uniform(const xmatrix<xcomplex<utype> >& a,       // real,imag in ([0,y],[0,y])
        const utype & y) {
      return uniform(a,utype(0.0),y,utype(0.0),y);
    }

  template<class utype> xmatrix<xcomplex<utype> >     // xcomplex xmatrix of uniforms
    uniform(const xmatrix<xcomplex<utype> >& a) {     // real,imag in ([0,1],[0,1])
      return uniform(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

// ----------------------------------------------------- gaussian distributions

namespace aurostd {
  // namespace aurostd
  template<class utype> xmatrix<utype>         // gaussian xmatrix<utype> (m,sigma)
    gaussian(const xmatrix<utype>& a,const utype& m,const utype& sigma)  {
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          //	a[i][j]=(utype) (utype) m+sqrt(sigma)*gaussian();
          a[i][j]=(utype) (utype) m+sigma*gaussian();
      return a;
    }

  template<class utype> xmatrix<utype>       // gaussian xmatrix<utype> (m=0,sigma)
    gaussian(const xmatrix<utype>& a,const utype& sigma)  {
      return gaussian(a,utype(0.0),utype(sigma));
    }

  template<class utype> xmatrix<utype>     // gaussian xmatrix<utype> (m=0,sigma=1)
    gaussian(const xmatrix<utype>& a)  {
      return gaussian(a,utype(0.0),utype(1.0));
    }

  template<class utype> xmatrix<xcomplex<utype> >    // xcomplex xmatrix of gaussians
    gaussian(const xmatrix<xcomplex<utype> >& a, //real,imag(m1,sigma1),(m2,sigma2)
        const utype & m1,const utype & sigma1,
        const utype & m2,const utype & sigma2) {
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          //	a[i][j]=xcomplex<utype>(m1+sqrt(sigma1)*gaussian(),
          //		       m2+sqrt(sigma2)*gaussian());
          a[i][j]=xcomplex<utype>(m1+sigma1*gaussian(),
              m2+sigma2*gaussian());
      return a;
    }

  template<class utype> xmatrix<xcomplex<utype> >    // xcomplex xmatrix of gaussians
    gaussian(const xmatrix<xcomplex<utype> >& a, // real,imag (0,sigma1),(0,sigma2)
        const utype & sigma1,const utype & sigma2) {
      return gaussian(a,utype(0.0),sigma1,utype(0.0),sigma2);
    }

  template<class utype> xmatrix<xcomplex<utype> >    // xcomplex xmatrix of gaussians
    gaussian(const xmatrix<xcomplex<utype> >& a,   // real,imag (0,sigma),(0,sigma)
        const utype & sigma) {
      return gaussian(a,utype(0.0),sigma,utype(0.0),sigma);
    }

  template<class utype> xmatrix<xcomplex<utype> >    // xcomplex xmatrix of gaussians
    gaussian(const xmatrix<xcomplex<utype> >& a) {      // real,imag in (0,1),(0,1)
      return gaussian(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

#endif

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
#ifdef __XTENSOR3_CPP

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor3<utype>              // uniform xtensor3<utype> [x,y]
    uniform(const xtensor3<utype>& a,const utype& x,const utype& y)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            a[i][j][k]=(utype) x+(y-x)*uniform(utype(1.0));
      return a;
    }

  template<class utype> xtensor3<utype>              // uniform xtensor3<utype> [0,y]
    uniform(const xtensor3<utype>& a,const utype& y)  {
      return uniform(a,utype(0.0),utype(y));
    }

  template<class utype> xtensor3<utype>              // uniform xtensor3<utype> [0,1]
    uniform(const xtensor3<utype>& a)  {
      return uniform(a,utype(0.0),utype(1.0));
    }


  template<class utype> xtensor3<xcomplex<utype> >     // xcomplex xtensor3 of uniforms
    uniform(const xtensor3<xcomplex<utype> >& a,   // real,imag in ([x1,y1],[x2,y2])
        const utype & x1,const utype & y1,
        const utype & x2,const utype & y2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            a[i][j][k]=xcomplex<utype>(x1+(y1-x1)*uniform(utype(1.0)),
                x2+(y2-x2)*uniform(utype(1.0)));
      return a;
    }

  template<class utype> xtensor3<xcomplex<utype> >     // xcomplex xtensor3 of uniforms
    uniform(const xtensor3<xcomplex<utype> >& a,       // real,imag in ([0,x],[0,y])
        const utype & x,const utype & y) {
      return uniform(a,utype(0.0),x,utype(0.0),y);
    }

  template<class utype> xtensor3<xcomplex<utype> >     // xcomplex xtensor3 of uniforms
    uniform(const xtensor3<xcomplex<utype> >& a,       // real,imag in ([0,y],[0,y])
        const utype & y) {
      return uniform(a,utype(0.0),y,utype(0.0),y);
    }

  template<class utype> xtensor3<xcomplex<utype> >     // xcomplex xtensor3 of uniforms
    uniform(const xtensor3<xcomplex<utype> >& a) {     // real,imag in ([0,1],[0,1])
      return uniform(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

// ----------------------------------------------------- gaussian distributions

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor3<utype>         // gaussian xtensor3<utype> (m,sigma)
    gaussian(const xtensor3<utype>& a,const utype& m,const utype& sigma)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            //	a[i][j][k]=(utype) (utype) m+sqrt(sigma)*gaussian();
            a[i][j][k]=(utype) (utype) m+sigma*gaussian();
      return a;
    }

  template<class utype> xtensor3<utype>       // gaussian xtensor3<utype> (m=0,sigma)
    gaussian(const xtensor3<utype>& a,const utype& sigma)  {
      return gaussian(a,utype(0.0),utype(sigma));
    }

  template<class utype> xtensor3<utype>     // gaussian xtensor3<utype> (m=0,sigma=1)
    gaussian(const xtensor3<utype>& a)  {
      return gaussian(a,utype(0.0),utype(1.0));
    }

  template<class utype> xtensor3<xcomplex<utype> >    // xcomplex xtensor3 of gaussians
    gaussian(const xtensor3<xcomplex<utype> >& a, //real,imag(m1,sigma1),(m2,sigma2)
        const utype & m1,const utype & sigma1,
        const utype & m2,const utype & sigma2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            //	a[i][j][k]=xcomplex<utype>(m1+sqrt(sigma1)*gaussian(),
            //		       m2+sqrt(sigma2)*gaussian());
            a[i][j][k]=xcomplex<utype>(m1+sigma1*gaussian(),
                m2+sigma2*gaussian());
      return a;
    }

  template<class utype> xtensor3<xcomplex<utype> >    // xcomplex xtensor3 of gaussians
    gaussian(const xtensor3<xcomplex<utype> >& a, // real,imag (0,sigma1),(0,sigma2)
        const utype & sigma1,const utype & sigma2) {
      return gaussian(a,utype(0.0),sigma1,utype(0.0),sigma2);
    }

  template<class utype> xtensor3<xcomplex<utype> >    // xcomplex xtensor3 of gaussians
    gaussian(const xtensor3<xcomplex<utype> >& a,   // real,imag (0,sigma),(0,sigma)
        const utype & sigma) {
      return gaussian(a,utype(0.0),sigma,utype(0.0),sigma);
    }

  template<class utype> xtensor3<xcomplex<utype> >    // xcomplex xtensor3 of gaussians
    gaussian(const xtensor3<xcomplex<utype> >& a) {      // real,imag in (0,1),(0,1)
      return gaussian(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

#endif
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
#ifdef __XTENSOR4_CPP

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor4<utype>              // uniform xtensor4<utype> [x,y]
    uniform(const xtensor4<utype>& a,const utype& x,const utype& y)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              a[i][j][k][l]=(utype) x+(y-x)*uniform(utype(1.0));
      return a;
    }

  template<class utype> xtensor4<utype>              // uniform xtensor4<utype> [0,y]
    uniform(const xtensor4<utype>& a,const utype& y)  {
      return uniform(a,utype(0.0),utype(y));
    }

  template<class utype> xtensor4<utype>              // uniform xtensor4<utype> [0,1]
    uniform(const xtensor4<utype>& a)  {
      return uniform(a,utype(0.0),utype(1.0));
    }


  template<class utype> xtensor4<xcomplex<utype> >     // xcomplex xtensor4 of uniforms
    uniform(const xtensor4<xcomplex<utype> >& a,   // real,imag in ([x1,y1],[x2,y2])
        const utype & x1,const utype & y1,
        const utype & x2,const utype & y2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              a[i][j][k][l]=xcomplex<utype>(x1+(y1-x1)*uniform(utype(1.0)),
                  x2+(y2-x2)*uniform(utype(1.0)));
      return a;
    }

  template<class utype> xtensor4<xcomplex<utype> >     // xcomplex xtensor4 of uniforms
    uniform(const xtensor4<xcomplex<utype> >& a,       // real,imag in ([0,x],[0,y])
        const utype & x,const utype & y) {
      return uniform(a,utype(0.0),x,utype(0.0),y);
    }

  template<class utype> xtensor4<xcomplex<utype> >     // xcomplex xtensor4 of uniforms
    uniform(const xtensor4<xcomplex<utype> >& a,       // real,imag in ([0,y],[0,y])
        const utype & y) {
      return uniform(a,utype(0.0),y,utype(0.0),y);
    }

  template<class utype> xtensor4<xcomplex<utype> >     // xcomplex xtensor4 of uniforms
    uniform(const xtensor4<xcomplex<utype> >& a) {     // real,imag in ([0,1],[0,1])
      return uniform(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

// ----------------------------------------------------- gaussian distributions

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor4<utype>         // gaussian xtensor4<utype> (m,sigma)
    gaussian(const xtensor4<utype>& a,const utype& m,const utype& sigma)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              //	a[i][j][k][l]=(utype) (utype) m+sqrt(sigma)*gaussian();
              a[i][j][k][l]=(utype) (utype) m+sigma*gaussian();
      return a;
    }

  template<class utype> xtensor4<utype>       // gaussian xtensor4<utype> (m=0,sigma)
    gaussian(const xtensor4<utype>& a,const utype& sigma)  {
      return gaussian(a,utype(0.0),utype(sigma));
    }

  template<class utype> xtensor4<utype>     // gaussian xtensor4<utype> (m=0,sigma=1)
    gaussian(const xtensor4<utype>& a)  {
      return gaussian(a,utype(0.0),utype(1.0));
    }

  template<class utype> xtensor4<xcomplex<utype> >    // xcomplex xtensor4 of gaussians
    gaussian(const xtensor4<xcomplex<utype> >& a, //real,imag(m1,sigma1),(m2,sigma2)
        const utype & m1,const utype & sigma1,
        const utype & m2,const utype & sigma2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              //	a[i][j][k][l]=xcomplex<utype>(m1+sqrt(sigma1)*gaussian(),
              //		       m2+sqrt(sigma2)*gaussian());
              \	    a[i][j][k][l]=xcomplex<utype>(m1+sigma1*gaussian(),
                  m2+sigma2*gaussian());
      return a;
    }

  template<class utype> xtensor4<xcomplex<utype> >    // xcomplex xtensor4 of gaussians
    gaussian(const xtensor4<xcomplex<utype> >& a, // real,imag (0,sigma1),(0,sigma2)
        const utype & sigma1,const utype & sigma2) {
      return gaussian(a,utype(0.0),sigma1,utype(0.0),sigma2);
    }

  template<class utype> xtensor4<xcomplex<utype> >    // xcomplex xtensor4 of gaussians
    gaussian(const xtensor4<xcomplex<utype> >& a,   // real,imag (0,sigma),(0,sigma)
        const utype & sigma) {
      return gaussian(a,utype(0.0),sigma,utype(0.0),sigma);
    }

  template<class utype> xtensor4<xcomplex<utype> >    // xcomplex xtensor4 of gaussians
    gaussian(const xtensor4<xcomplex<utype> >& a) {      // real,imag in (0,1),(0,1)
      return gaussian(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

#endif
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
#ifdef __XTENSOR5_CPP

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor5<utype>              // uniform xtensor5<utype> [x,y]
    uniform(const xtensor5<utype>& a,const utype& x,const utype& y)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                a[i][j][k][l][m]=(utype) x+(y-x)*uniform(utype(1.0));
      return a;
    }

  template<class utype> xtensor5<utype>              // uniform xtensor5<utype> [0,y]
    uniform(const xtensor5<utype>& a,const utype& y)  {
      return uniform(a,utype(0.0),utype(y));
    }

  template<class utype> xtensor5<utype>              // uniform xtensor5<utype> [0,1]
    uniform(const xtensor5<utype>& a)  {
      return uniform(a,utype(0.0),utype(1.0));
    }


  template<class utype> xtensor5<xcomplex<utype> >     // xcomplex xtensor5 of uniforms
    uniform(const xtensor5<xcomplex<utype> >& a,   // real,imag in ([x1,y1],[x2,y2])
        const utype & x1,const utype & y1,
        const utype & x2,const utype & y2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                a[i][j][k][l][m]=xcomplex<utype>(x1+(y1-x1)*uniform(utype(1.0)),
                    x2+(y2-x2)*uniform(utype(1.0)));
      return a;
    }

  template<class utype> xtensor5<xcomplex<utype> >     // xcomplex xtensor5 of uniforms
    uniform(const xtensor5<xcomplex<utype> >& a,       // real,imag in ([0,x],[0,y])
        const utype & x,const utype & y) {
      return uniform(a,utype(0.0),x,utype(0.0),y);
    }

  template<class utype> xtensor5<xcomplex<utype> >     // xcomplex xtensor5 of uniforms
    uniform(const xtensor5<xcomplex<utype> >& a,       // real,imag in ([0,y],[0,y])
        const utype & y) {
      return uniform(a,utype(0.0),y,utype(0.0),y);
    }

  template<class utype> xtensor5<xcomplex<utype> >     // xcomplex xtensor5 of uniforms
    uniform(const xtensor5<xcomplex<utype> >& a) {     // real,imag in ([0,1],[0,1])
      return uniform(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

// ----------------------------------------------------- gaussian distributions

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor5<utype>         // gaussian xtensor5<utype> (m,sigma)
    gaussian(const xtensor5<utype>& a,const utype& m,const utype& sigma)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                //	a[i][j][k][l][m]=(utype) (utype) m+sqrt(sigma)*gaussian();
                a[i][j][k][l][m]=(utype) (utype) m+sigma*gaussian();
      return a;
    }

  template<class utype> xtensor5<utype>       // gaussian xtensor5<utype> (m=0,sigma)
    gaussian(const xtensor5<utype>& a,const utype& sigma)  {
      return gaussian(a,utype(0.0),utype(sigma));
    }

  template<class utype> xtensor5<utype>     // gaussian xtensor5<utype> (m=0,sigma=1)
    gaussian(const xtensor5<utype>& a)  {
      return gaussian(a,utype(0.0),utype(1.0));
    }

  template<class utype> xtensor5<xcomplex<utype> >    // xcomplex xtensor5 of gaussians
    gaussian(const xtensor5<xcomplex<utype> >& a, //real,imag(m1,sigma1),(m2,sigma2)
        const utype & m1,const utype & sigma1,
        const utype & m2,const utype & sigma2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                //	a[i][j][k][l][m]=xcomplex<utype>(m1+sqrt(sigma1)*gaussian(),
                //		       m2+sqrt(sigma2)*gaussian());
                a[i][j][k][l][m]=xcomplex<utype>(m1+sigma1*gaussian(),
                    m2+sigma2*gaussian());
      return a;
    }

  template<class utype> xtensor5<xcomplex<utype> >    // xcomplex xtensor5 of gaussians
    gaussian(const xtensor5<xcomplex<utype> >& a, // real,imag (0,sigma1),(0,sigma2)
        const utype & sigma1,const utype & sigma2) {
      return gaussian(a,utype(0.0),sigma1,utype(0.0),sigma2);
    }

  template<class utype> xtensor5<xcomplex<utype> >    // xcomplex xtensor5 of gaussians
    gaussian(const xtensor5<xcomplex<utype> >& a,   // real,imag (0,sigma),(0,sigma)
        const utype & sigma) {
      return gaussian(a,utype(0.0),sigma,utype(0.0),sigma);
    }

  template<class utype> xtensor5<xcomplex<utype> >    // xcomplex xtensor5 of gaussians
    gaussian(const xtensor5<xcomplex<utype> >& a) {      // real,imag in (0,1),(0,1)
      return gaussian(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

#endif
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
#ifdef __XTENSOR6_CPP

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor6<utype>              // uniform xtensor6<utype> [x,y]
    uniform(const xtensor6<utype>& a,const utype& x,const utype& y)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  a[i][j][k][l][m][n]=(utype) x+(y-x)*uniform(utype(1.0));
      return a;
    }

  template<class utype> xtensor6<utype>              // uniform xtensor6<utype> [0,y]
    uniform(const xtensor6<utype>& a,const utype& y)  {
      return uniform(a,utype(0.0),utype(y));
    }

  template<class utype> xtensor6<utype>              // uniform xtensor6<utype> [0,1]
    uniform(const xtensor6<utype>& a)  {
      return uniform(a,utype(0.0),utype(1.0));
    }


  template<class utype> xtensor6<xcomplex<utype> >     // xcomplex xtensor6 of uniforms
    uniform(const xtensor6<xcomplex<utype> >& a,   // real,imag in ([x1,y1],[x2,y2])
        const utype & x1,const utype & y1,
        const utype & x2,const utype & y2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  a[i][j][k][l][m][n]=xcomplex<utype>(x1+(y1-x1)*uniform(utype(1.0)),
                      x2+(y2-x2)*uniform(utype(1.0)));
      return a;
    }

  template<class utype> xtensor6<xcomplex<utype> >     // xcomplex xtensor6 of uniforms
    uniform(const xtensor6<xcomplex<utype> >& a,       // real,imag in ([0,x],[0,y])
        const utype & x,const utype & y) {
      return uniform(a,utype(0.0),x,utype(0.0),y);
    }

  template<class utype> xtensor6<xcomplex<utype> >     // xcomplex xtensor6 of uniforms
    uniform(const xtensor6<xcomplex<utype> >& a,       // real,imag in ([0,y],[0,y])
        const utype & y) {
      return uniform(a,utype(0.0),y,utype(0.0),y);
    }

  template<class utype> xtensor6<xcomplex<utype> >     // xcomplex xtensor6 of uniforms
    uniform(const xtensor6<xcomplex<utype> >& a) {     // real,imag in ([0,1],[0,1])
      return uniform(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

// ----------------------------------------------------- gaussian distributions

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor6<utype>         // gaussian xtensor6<utype> (m,sigma)
    gaussian(const xtensor6<utype>& a,const utype& m,const utype& sigma)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  //	a[i][j][k][l][m][n]=(utype) (utype) m+sqrt(sigma)*gaussian();
                  a[i][j][k][l][m][n]=(utype) (utype) m+sigma*gaussian();
      return a;
    }

  template<class utype> xtensor6<utype>       // gaussian xtensor6<utype> (m=0,sigma)
    gaussian(const xtensor6<utype>& a,const utype& sigma)  {
      return gaussian(a,utype(0.0),utype(sigma));
    }

  template<class utype> xtensor6<utype>     // gaussian xtensor6<utype> (m=0,sigma=1)
    gaussian(const xtensor6<utype>& a)  {
      return gaussian(a,utype(0.0),utype(1.0));
    }

  template<class utype> xtensor6<xcomplex<utype> >    // xcomplex xtensor6 of gaussians
    gaussian(const xtensor6<xcomplex<utype> >& a, //real,imag(m1,sigma1),(m2,sigma2)
        const utype & m1,const utype & sigma1,
        const utype & m2,const utype & sigma2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  //	a[i][j][k][l][m][n]=xcomplex<utype>(m1+sqrt(sigma1)*gaussian(),
                  //		       m2+sqrt(sigma2)*gaussian());
                  a[i][j][k][l][m][n]=xcomplex<utype>(m1+sigma1*gaussian(),
                      m2+sigma2*gaussian());
      return a;
    }

  template<class utype> xtensor6<xcomplex<utype> >    // xcomplex xtensor6 of gaussians
    gaussian(const xtensor6<xcomplex<utype> >& a, // real,imag (0,sigma1),(0,sigma2)
        const utype & sigma1,const utype & sigma2) {
      return gaussian(a,utype(0.0),sigma1,utype(0.0),sigma2);
    }

  template<class utype> xtensor6<xcomplex<utype> >    // xcomplex xtensor6 of gaussians
    gaussian(const xtensor6<xcomplex<utype> >& a,   // real,imag (0,sigma),(0,sigma)
        const utype & sigma) {
      return gaussian(a,utype(0.0),sigma,utype(0.0),sigma);
    }

  template<class utype> xtensor6<xcomplex<utype> >    // xcomplex xtensor6 of gaussians
    gaussian(const xtensor6<xcomplex<utype> >& a) {      // real,imag in (0,1),(0,1)
      return gaussian(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

#endif
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
#ifdef __XTENSOR7_CPP

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor7<utype>              // uniform xtensor7<utype> [x,y]
    uniform(const xtensor7<utype>& a,const utype& x,const utype& y)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  for(int o=a.lindex[7];o<=a.uindex[7];o++)
                    a[i][j][k][l][m][n][o]=(utype) x+(y-x)*uniform(utype(1.0));
      return a;
    }

  template<class utype> xtensor7<utype>              // uniform xtensor7<utype> [0,y]
    uniform(const xtensor7<utype>& a,const utype& y)  {
      return uniform(a,utype(0.0),utype(y));
    }

  template<class utype> xtensor7<utype>              // uniform xtensor7<utype> [0,1]
    uniform(const xtensor7<utype>& a)  {
      return uniform(a,utype(0.0),utype(1.0));
    }


  template<class utype> xtensor7<xcomplex<utype> >     // xcomplex xtensor7 of uniforms
    uniform(const xtensor7<xcomplex<utype> >& a,   // real,imag in ([x1,y1],[x2,y2])
        const utype & x1,const utype & y1,
        const utype & x2,const utype & y2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  for(int o=a.lindex[7];o<=a.uindex[7];o++)
                    a[i][j][k][l][m][n][o]=xcomplex<utype>(x1+(y1-x1)*uniform(utype(1.0)),
                        x2+(y2-x2)*uniform(utype(1.0)));
      return a;
    }

  template<class utype> xtensor7<xcomplex<utype> >     // xcomplex xtensor7 of uniforms
    uniform(const xtensor7<xcomplex<utype> >& a,       // real,imag in ([0,x],[0,y])
        const utype & x,const utype & y) {
      return uniform(a,utype(0.0),x,utype(0.0),y);
    }

  template<class utype> xtensor7<xcomplex<utype> >     // xcomplex xtensor7 of uniforms
    uniform(const xtensor7<xcomplex<utype> >& a,       // real,imag in ([0,y],[0,y])
        const utype & y) {
      return uniform(a,utype(0.0),y,utype(0.0),y);
    }

  template<class utype> xtensor7<xcomplex<utype> >     // xcomplex xtensor7 of uniforms
    uniform(const xtensor7<xcomplex<utype> >& a) {     // real,imag in ([0,1],[0,1])
      return uniform(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

// ----------------------------------------------------- gaussian distributions

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor7<utype>         // gaussian xtensor7<utype> (m,sigma)
    gaussian(const xtensor7<utype>& a,const utype& m,const utype& sigma)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  for(int o=a.lindex[7];o<=a.uindex[7];o++)
                    //	a[i][j][k][l][m][n][o]=(utype) (utype) m+sqrt(sigma)*gaussian();
                    a[i][j][k][l][m][n][o]=(utype) (utype) m+sigma*gaussian();
      return a;
    }

  template<class utype> xtensor7<utype>       // gaussian xtensor7<utype> (m=0,sigma)
    gaussian(const xtensor7<utype>& a,const utype& sigma)  {
      return gaussian(a,utype(0.0),utype(sigma));
    }

  template<class utype> xtensor7<utype>     // gaussian xtensor7<utype> (m=0,sigma=1)
    gaussian(const xtensor7<utype>& a)  {
      return gaussian(a,utype(0.0),utype(1.0));
    }

  template<class utype> xtensor7<xcomplex<utype> >    // xcomplex xtensor7 of gaussians
    gaussian(const xtensor7<xcomplex<utype> >& a, //real,imag(m1,sigma1),(m2,sigma2)
        const utype & m1,const utype & sigma1,
        const utype & m2,const utype & sigma2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  for(int o=a.lindex[7];o<=a.uindex[7];o++)
                    //	a[i][j][k][l][m][n][o]=xcomplex<utype>(m1+sqrt(sigma1)*gaussian(),
                    //		       m2+sqrt(sigma2)*gaussian());
                    a[i][j][k][l][m][n][o]=xcomplex<utype>(m1+sigma1*gaussian(),
                        m2+sigma2*gaussian());
      return a;
    }

  template<class utype> xtensor7<xcomplex<utype> >    // xcomplex xtensor7 of gaussians
    gaussian(const xtensor7<xcomplex<utype> >& a, // real,imag (0,sigma1),(0,sigma2)
        const utype & sigma1,const utype & sigma2) {
      return gaussian(a,utype(0.0),sigma1,utype(0.0),sigma2);
    }

  template<class utype> xtensor7<xcomplex<utype> >    // xcomplex xtensor7 of gaussians
    gaussian(const xtensor7<xcomplex<utype> >& a,   // real,imag (0,sigma),(0,sigma)
        const utype & sigma) {
      return gaussian(a,utype(0.0),sigma,utype(0.0),sigma);
    }

  template<class utype> xtensor7<xcomplex<utype> >    // xcomplex xtensor7 of gaussians
    gaussian(const xtensor7<xcomplex<utype> >& a) {      // real,imag in (0,1),(0,1)
      return gaussian(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

#endif
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
#ifdef __XTENSOR8_CPP

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor8<utype>              // uniform xtensor8<utype> [x,y]
    uniform(const xtensor8<utype>& a,const utype& x,const utype& y)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  for(int o=a.lindex[7];o<=a.uindex[7];o++)
                    for(int p=a.lindex[8];p<=a.uindex[8];p++)
                      a[i][j][k][l][m][n][o][p]=(utype) x+(y-x)*uniform(utype(1.0));
      return a;
    }

  template<class utype> xtensor8<utype>              // uniform xtensor8<utype> [0,y]
    uniform(const xtensor8<utype>& a,const utype& y)  {
      return uniform(a,utype(0.0),utype(y));
    }

  template<class utype> xtensor8<utype>              // uniform xtensor8<utype> [0,1]
    uniform(const xtensor8<utype>& a)  {
      return uniform(a,utype(0.0),utype(1.0));
    }


  template<class utype> xtensor8<xcomplex<utype> >     // xcomplex xtensor8 of uniforms
    uniform(const xtensor8<xcomplex<utype> >& a,   // real,imag in ([x1,y1],[x2,y2])
        const utype & x1,const utype & y1,
        const utype & x2,const utype & y2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  for(int o=a.lindex[7];o<=a.uindex[7];o++)
                    for(int p=a.lindex[8];p<=a.uindex[8];p++)
                      a[i][j][k][l][m][n][o][p]=xcomplex<utype>(x1+(y1-x1)*uniform(utype(1.0)),
                          x2+(y2-x2)*uniform(utype(1.0)));
      return a;
    }

  template<class utype> xtensor8<xcomplex<utype> >     // xcomplex xtensor8 of uniforms
    uniform(const xtensor8<xcomplex<utype> >& a,       // real,imag in ([0,x],[0,y])
        const utype & x,const utype & y) {
      return uniform(a,utype(0.0),x,utype(0.0),y);
    }

  template<class utype> xtensor8<xcomplex<utype> >     // xcomplex xtensor8 of uniforms
    uniform(const xtensor8<xcomplex<utype> >& a,       // real,imag in ([0,y],[0,y])
        const utype & y) {
      return uniform(a,utype(0.0),y,utype(0.0),y);
    }

  template<class utype> xtensor8<xcomplex<utype> >     // xcomplex xtensor8 of uniforms
    uniform(const xtensor8<xcomplex<utype> >& a) {     // real,imag in ([0,1],[0,1])
      return uniform(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

// ----------------------------------------------------- gaussian distributions

namespace aurostd {
  // namespace aurostd
  template<class utype> xtensor8<utype>         // gaussian xtensor8<utype> (m,sigma)
    gaussian(const xtensor8<utype>& a,const utype& m,const utype& sigma)  {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  for(int o=a.lindex[7];o<=a.uindex[7];o++)
                    for(int p=a.lindex[8];p<=a.uindex[8];p++)
                      //	a[i][j][k][l][m][n][o][p]=(utype) (utype) m+sqrt(sigma)*gaussian();
                      a[i][j][k][l][m][n][o][p]=(utype) (utype) m+sigma*gaussian();
      return a;
    }

  template<class utype> xtensor8<utype>       // gaussian xtensor8<utype> (m=0,sigma)
    gaussian(const xtensor8<utype>& a,const utype& sigma)  {
      return gaussian(a,utype(0.0),utype(sigma));
    }

  template<class utype> xtensor8<utype>     // gaussian xtensor8<utype> (m=0,sigma=1)
    gaussian(const xtensor8<utype>& a)  {
      return gaussian(a,utype(0.0),utype(1.0));
    }

  template<class utype> xtensor8<xcomplex<utype> >    // xcomplex xtensor8 of gaussians
    gaussian(const xtensor8<xcomplex<utype> >& a, //real,imag(m1,sigma1),(m2,sigma2)
        const utype & m1,const utype & sigma1,
        const utype & m2,const utype & sigma2) {
      for(int i=a.lindex[1];i<=a.uindex[1];i++)
        for(int j=a.lindex[2];j<=a.uindex[2];j++)
          for(int k=a.lindex[3];k<=a.uindex[3];k++)
            for(int l=a.lindex[4];l<=a.uindex[4];l++)
              for(int m=a.lindex[5];m<=a.uindex[5];m++)
                for(int n=a.lindex[6];n<=a.uindex[6];n++)
                  for(int o=a.lindex[7];o<=a.uindex[7];o++)
                    for(int p=a.lindex[8];p<=a.uindex[8];p++)
                      //	a[i][j][k][l][m][n][o][p]=xcomplex<utype>(m1+sqrt(sigma1)*gaussian(),
                      //		       m2+sqrt(sigma2)*gaussian());
                      a[i][j][k][l][m][n][o][p]=xcomplex<utype>(m1+sigma1*gaussian(),
                          m2+sigma2*gaussian());
      return a;
    }

  template<class utype> xtensor8<xcomplex<utype> >    // xcomplex xtensor8 of gaussians
    gaussian(const xtensor8<xcomplex<utype> >& a, // real,imag (0,sigma1),(0,sigma2)
        const utype & sigma1,const utype & sigma2) {
      return gaussian(a,utype(0.0),sigma1,utype(0.0),sigma2);
    }

  template<class utype> xtensor8<xcomplex<utype> >    // xcomplex xtensor8 of gaussians
    gaussian(const xtensor8<xcomplex<utype> >& a,   // real,imag (0,sigma),(0,sigma)
        const utype & sigma) {
      return gaussian(a,utype(0.0),sigma,utype(0.0),sigma);
    }

  template<class utype> xtensor8<xcomplex<utype> >    // xcomplex xtensor8 of gaussians
    gaussian(const xtensor8<xcomplex<utype> >& a) {      // real,imag in (0,1),(0,1)
      return gaussian(a,utype(0.0),utype(1.0),utype(0.0),utype(1.0));
    }
}

#endif

// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XRANDOM_IMPLEMENTATIONS_



// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************

