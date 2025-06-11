// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XRANDOM_H_
#define _AUROSTD_XRANDOM_H_

#include "aurostd_defs.h"

namespace aurostd {
  // namespace aurostd
  // ----------------------------------------------------------------------------
  // ---------------------------------- random function and procedures on numbers

  long int _random_initialize(long int _num);               // initialize randoms
  long int _random_initialize();// initialize randoms with clock based number
  // ---------------------------------------------------------- uniform deviation
  double ran0();                                          // uniform in [0,1]
  double ran1();                                          // uniform in [0,1]
  double ran2();                                          // uniform in [0,1]
  template <class utype>
  utype uniform(const utype &)    // uniform in [0,x]
      __xprototype;
  template <class utype>
  inline utype uniform(const utype &,    // uniform in [x,y]
                       const utype &) __xprototype;
  // --------------------------------------------------------- gaussian deviation
  double gaussian();                                     // gaussian(m=0,v=1)
  template <class utype>
  inline utype gaussian(const utype &)    // gaussian(m=0,v)
      __xprototype;
  template <class utype>
  inline utype gaussian(const utype &,      // gaussian(m,v)
                        const utype &) __xprototype;
  // ------------------------------------------------------ exponential deviation
  inline double expdev()                              // exponential lambda=1
    // returns a exponentially distribuited deviate with lambda=1
      __xprototype;
  template <class utype>
  inline utype expdev(const utype &)   // exponential lambda
      __xprototype;

  // ---------------------------------------------------------- laplace deviation
  inline double laplacedev()                            // laplace mean=0 a=1
    // returns a laplace distribuited deviate p(x)=1/2a*exp(-|x-lambda|/a)
    // sigma^2=2a*a mean=lambda
      __xprototype;

  template <class utype>
  inline utype                          // laplace mean=0 a
  laplacedev(const utype &)
    // returns a laplace distribuited deviate p(x)=1/2a*exp(-|x-lambda|/a)
      __xprototype;

  template <class utype>
  inline utype                            // laplace mean a
  laplacedev(const utype &a, const utype &lambda)
    // returns a laplace distribuited deviate p(x)=1/2a*exp(-|x-lambda|/a)
      __xprototype;
} // namespace aurostd

// ---------------------------------------------------------- poisson deviation

// ????????? not yet done

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// -------------------------- uniform function and procedures on xvector, xmatrix

#ifdef _AUROSTD_XVECTOR_CPP_

#include "aurostd_xvector.h"

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xvector<utype>                 // xvector of uniforms [0,1]
  uniform(const xvector<utype> &) __xprototype;

  template <class utype>
  xvector<utype>                 // xvector of uniforms [0,y]
  uniform(const xvector<utype> &, const utype &) __xprototype;

  template <class utype>
  xvector<utype>                 // xvector of uniforms [x,y]
  uniform(const xvector<utype> &, const utype &, const utype &) __xprototype;

  template <class utype>
  xvector<std::complex<utype>>     // complex xvector of uniforms
  uniform(const xvector<std::complex<utype>> &)         // real,imag in ([0,1],[0,1])
      __xprototype;

  template <class utype>
  xvector<std::complex<utype>>     // complex xvector of uniforms
  uniform(const xvector<std::complex<utype>> &,         // real,imag in ([0,y],[0,y])
          const utype &) __xprototype;

  template <class utype>
  xvector<std::complex<utype>>     // complex xvector of uniforms
  uniform(const xvector<std::complex<utype>> &,       // real,imag in ([0,y1],[0,y2])
          const utype &,
          const utype &) __xprototype;

  template <class utype>
  xvector<std::complex<utype>>     // complex xvector of uniforms
  uniform(const xvector<std::complex<utype>> &,     // real,imag in ([x1,y1],[x2,y2])
          const utype &,
          const utype &,
          const utype &,
          const utype &) __xprototype;
} // namespace aurostd
#endif

// ---------------------------------- uniform function and procedures on xmatrix

#ifdef _AUROSTD_XMATRIX_CPP_

#include "aurostd_xmatrix.h"

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xmatrix<utype>                 // xmatrix of uniforms [0,1]
  uniform(const xmatrix<utype> &) __xprototype;

  template <class utype>
  xmatrix<utype>                 // xmatrix of uniforms [0,y]
  uniform(const xmatrix<utype> &, const utype &) __xprototype;

  template <class utype>
  xmatrix<utype>                 // xmatrix of uniforms [x,y]
  uniform(const xmatrix<utype> &, const utype &, const utype &) __xprototype;

  template <class utype>
  xmatrix<std::complex<utype>>     // complex xmatrix of uniforms
  uniform(const xmatrix<std::complex<utype>> &)         // real,imag in ([0,1],[0,1])
      __xprototype;

  template <class utype>
  xmatrix<std::complex<utype>>     // complex xmatrix of uniforms
  uniform(const xmatrix<std::complex<utype>> &,         // real,imag in ([0,y],[0,y])
          const utype &) __xprototype;

  template <class utype>
  xmatrix<std::complex<utype>>     // complex xmatrix of uniforms
  uniform(const xmatrix<std::complex<utype>> &,         // real,imag in ([0,x],[0,y])
          const utype &,
          const utype &) __xprototype;

  template <class utype>
  xmatrix<std::complex<utype>>     // complex xmatrix of uniforms
  uniform(const xmatrix<std::complex<utype>> &,     // real,imag in ([x1,y1],[x2,y2])
          const utype &,
          const utype &,
          const utype &,
          const utype &) __xprototype;
} // namespace aurostd

// --------------------------------- gaussian function and procedures on xmatrix

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xmatrix<utype>             // gaussian xmatrix<utype> (m,v)
  gaussian(const xmatrix<utype> &a, const utype &m, const utype &v) __xprototype;

  template <class utype>
  xmatrix<utype>           // gaussian xmatrix<utype> (m=0,v)
  gaussian(const xmatrix<utype> &a, const utype &v) __xprototype;

  template <class utype>
  xmatrix<utype>         // gaussian xmatrix<utype> (m=0,v=1)
  gaussian(const xmatrix<utype> &a) __xprototype;

  template <class utype>
  xmatrix<std::complex<utype>>    // complex xmatrix of gaussians
  gaussian(const xmatrix<std::complex<utype>> &a,    // real,imag in (m1,v1),(m2,v2)
           const utype &m2,
           const utype &v2) __xprototype;

  template <class utype>
  xmatrix<std::complex<utype>>    // complex xmatrix of gaussians
  gaussian(const xmatrix<std::complex<utype>> &a,      // real,imag in (0,v1),(0,v2)
           const utype &v1,
           const utype &v2) __xprototype;

  template <class utype>
  xmatrix<std::complex<utype>>    // complex xmatrix of gaussians
  gaussian(const xmatrix<std::complex<utype>> &a,        // real,imag in (0,v),(0,v)
           const utype &v) __xprototype;

  template <class utype>
  xmatrix<std::complex<utype>>    // complex xmatrix of gaussians
  gaussian(const xmatrix<std::complex<utype>> &a)       // real,imag in (0,1),(0,1)
      __xprototype;
} // namespace aurostd
#endif

// ---------------------------------- uniform function and procedures on xtensor3

#ifdef __XTENSOR3_CPP

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor3<utype>                 // xtensor3 of uniforms [0,1]
  uniform(const xtensor3<utype> &) __xprototype;

  template <class utype>
  xtensor3<utype>                 // xtensor3 of uniforms [0,y]
  uniform(const xtensor3<utype> &, const utype &) __xprototype;

  template <class utype>
  xtensor3<utype>                 // xtensor3 of uniforms [x,y]
  uniform(const xtensor3<utype> &, const utype &, const utype &) __xprototype;

  template <class utype>
  xtensor3<std::complex<utype>>     // complex xtensor3 of uniforms
  uniform(const xtensor3<std::complex<utype>> &)         // real,imag in ([0,1],[0,1])
      __xprototype;

  template <class utype>
  xtensor3<std::complex<utype>>     // complex xtensor3 of uniforms
  uniform(const xtensor3<std::complex<utype>> &,         // real,imag in ([0,y],[0,y])
          const utype &) __xprototype;

  template <class utype>
  xtensor3<std::complex<utype>>     // complex xtensor3 of uniforms
  uniform(const xtensor3<std::complex<utype>> &,         // real,imag in ([0,x],[0,y])
          const utype &,
          const utype &) __xprototype;

  template <class utype>
  xtensor3<std::complex<utype>>     // complex xtensor3 of uniforms
  uniform(const xtensor3<std::complex<utype>> &,     // real,imag in ([x1,y1],[x2,y2])
          const utype &,
          const utype &,
          const utype &,
          const utype &) __xprototype;
} // namespace aurostd

// --------------------------------- gaussian function and procedures on xtensor3

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor3<utype>             // gaussian xtensor3<utype> (m,v)
  gaussian(const xtensor3<utype> &a, const utype &m, const utype &v) __xprototype;

  template <class utype>
  xtensor3<utype>           // gaussian xtensor3<utype> (m=0,v)
  gaussian(const xtensor3<utype> &a, const utype &v) __xprototype;

  template <class utype>
  xtensor3<utype>         // gaussian xtensor3<utype> (m=0,v=1)
  gaussian(const xtensor3<utype> &a) __xprototype;

  template <class utype>
  xtensor3<std::complex<utype>>    // complex xtensor3 of gaussians
  gaussian(const xtensor3<std::complex<utype>> &a,    // real,imag in (m1,v1),(m2,v2)
           const utype &m2,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor3<std::complex<utype>>    // complex xtensor3 of gaussians
  gaussian(const xtensor3<std::complex<utype>> &a,      // real,imag in (0,v1),(0,v2)
           const utype &v1,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor3<std::complex<utype>>    // complex xtensor3 of gaussians
  gaussian(const xtensor3<std::complex<utype>> &a,        // real,imag in (0,v),(0,v)
           const utype &v) __xprototype;

  template <class utype>
  xtensor3<std::complex<utype>>    // complex xtensor3 of gaussians
  gaussian(const xtensor3<std::complex<utype>> &a)       // real,imag in (0,1),(0,1)
      __xprototype;
} // namespace aurostd

#endif //__XTENSOR3_CPP

// ---------------------------------- uniform function and procedures on xtensor4

#ifdef __XTENSOR4_CPP

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor4<utype>                 // xtensor4 of uniforms [0,1]
  uniform(const xtensor4<utype> &) __xprototype;

  template <class utype>
  xtensor4<utype>                 // xtensor4 of uniforms [0,y]
  uniform(const xtensor4<utype> &, const utype &) __xprototype;

  template <class utype>
  xtensor4<utype>                 // xtensor4 of uniforms [x,y]
  uniform(const xtensor4<utype> &, const utype &, const utype &) __xprototype;

  template <class utype>
  xtensor4<std::complex<utype>>     // complex xtensor4 of uniforms
  uniform(const xtensor4<std::complex<utype>> &)         // real,imag in ([0,1],[0,1])
      __xprototype;

  template <class utype>
  xtensor4<std::complex<utype>>     // complex xtensor4 of uniforms
  uniform(const xtensor4<std::complex<utype>> &,         // real,imag in ([0,y],[0,y])
          const utype &) __xprototype;

  template <class utype>
  xtensor4<std::complex<utype>>     // complex xtensor4 of uniforms
  uniform(const xtensor4<std::complex<utype>> &,         // real,imag in ([0,x],[0,y])
          const utype &,
          const utype &) __xprototype;

  template <class utype>
  xtensor4<std::complex<utype>>     // complex xtensor4 of uniforms
  uniform(const xtensor4<std::complex<utype>> &,     // real,imag in ([x1,y1],[x2,y2])
          const utype &,
          const utype &,
          const utype &,
          const utype &) __xprototype;
} // namespace aurostd

// --------------------------------- gaussian function and procedures on xtensor4

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor4<utype>             // gaussian xtensor4<utype> (m,v)
  gaussian(const xtensor4<utype> &a, const utype &m, const utype &v) __xprototype;

  template <class utype>
  xtensor4<utype>           // gaussian xtensor4<utype> (m=0,v)
  gaussian(const xtensor4<utype> &a, const utype &v) __xprototype;

  template <class utype>
  xtensor4<utype>         // gaussian xtensor4<utype> (m=0,v=1)
  gaussian(const xtensor4<utype> &a) __xprototype;

  template <class utype>
  xtensor4<std::complex<utype>>    // complex xtensor4 of gaussians
  gaussian(const xtensor4<std::complex<utype>> &a,    // real,imag in (m1,v1),(m2,v2)
           const utype &m2,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor4<std::complex<utype>>    // complex xtensor4 of gaussians
  gaussian(const xtensor4<std::complex<utype>> &a,      // real,imag in (0,v1),(0,v2)
           const utype &v1,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor4<std::complex<utype>>    // complex xtensor4 of gaussians
  gaussian(const xtensor4<std::complex<utype>> &a,        // real,imag in (0,v),(0,v)
           const utype &v) __xprototype;

  template <class utype>
  xtensor4<std::complex<utype>>    // complex xtensor4 of gaussians
  gaussian(const xtensor4<std::complex<utype>> &a)       // real,imag in (0,1),(0,1)
      __xprototype;
} // namespace aurostd

#endif //__XTENSOR4_CPP

// ---------------------------------- uniform function and procedures on xtensor5

#ifdef __XTENSOR5_CPP

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor5<utype>                 // xtensor5 of uniforms [0,1]
  uniform(const xtensor5<utype> &) __xprototype;

  template <class utype>
  xtensor5<utype>                 // xtensor5 of uniforms [0,y]
  uniform(const xtensor5<utype> &, const utype &) __xprototype;

  template <class utype>
  xtensor5<utype>                 // xtensor5 of uniforms [x,y]
  uniform(const xtensor5<utype> &, const utype &, const utype &) __xprototype;

  template <class utype>
  xtensor5<std::complex<utype>>     // complex xtensor5 of uniforms
  uniform(const xtensor5<std::complex<utype>> &)         // real,imag in ([0,1],[0,1])
      __xprototype;

  template <class utype>
  xtensor5<std::complex<utype>>     // complex xtensor5 of uniforms
  uniform(const xtensor5<std::complex<utype>> &,         // real,imag in ([0,y],[0,y])
          const utype &) __xprototype;

  template <class utype>
  xtensor5<std::complex<utype>>     // complex xtensor5 of uniforms
  uniform(const xtensor5<std::complex<utype>> &,         // real,imag in ([0,x],[0,y])
          const utype &,
          const utype &) __xprototype;

  template <class utype>
  xtensor5<std::complex<utype>>     // complex xtensor5 of uniforms
  uniform(const xtensor5<std::complex<utype>> &,     // real,imag in ([x1,y1],[x2,y2])
          const utype &,
          const utype &,
          const utype &,
          const utype &) __xprototype;
} // namespace aurostd

// --------------------------------- gaussian function and procedures on xtensor5

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor5<utype>             // gaussian xtensor5<utype> (m,v)
  gaussian(const xtensor5<utype> &a, const utype &m, const utype &v) __xprototype;

  template <class utype>
  xtensor5<utype>           // gaussian xtensor5<utype> (m=0,v)
  gaussian(const xtensor5<utype> &a, const utype &v) __xprototype;

  template <class utype>
  xtensor5<utype>         // gaussian xtensor5<utype> (m=0,v=1)
  gaussian(const xtensor5<utype> &a) __xprototype;

  template <class utype>
  xtensor5<std::complex<utype>>    // complex xtensor5 of gaussians
  gaussian(const xtensor5<std::complex<utype>> &a,    // real,imag in (m1,v1),(m2,v2)
           const utype &m2,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor5<std::complex<utype>>    // complex xtensor5 of gaussians
  gaussian(const xtensor5<std::complex<utype>> &a,      // real,imag in (0,v1),(0,v2)
           const utype &v1,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor5<std::complex<utype>>    // complex xtensor5 of gaussians
  gaussian(const xtensor5<std::complex<utype>> &a,        // real,imag in (0,v),(0,v)
           const utype &v) __xprototype;

  template <class utype>
  xtensor5<std::complex<utype>>    // complex xtensor5 of gaussians
  gaussian(const xtensor5<std::complex<utype>> &a)       // real,imag in (0,1),(0,1)
      __xprototype;
} // namespace aurostd

#endif //__XTENSOR5_CPP

// ---------------------------------- uniform function and procedures on xtensor6

#ifdef __XTENSOR6_CPP

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor6<utype>                 // xtensor6 of uniforms [0,1]
  uniform(const xtensor6<utype> &) __xprototype;

  template <class utype>
  xtensor6<utype>                 // xtensor6 of uniforms [0,y]
  uniform(const xtensor6<utype> &, const utype &) __xprototype;

  template <class utype>
  xtensor6<utype>                 // xtensor6 of uniforms [x,y]
  uniform(const xtensor6<utype> &, const utype &, const utype &) __xprototype;

  template <class utype>
  xtensor6<std::complex<utype>>     // complex xtensor6 of uniforms
  uniform(const xtensor6<std::complex<utype>> &)         // real,imag in ([0,1],[0,1])
      __xprototype;

  template <class utype>
  xtensor6<std::complex<utype>>     // complex xtensor6 of uniforms
  uniform(const xtensor6<std::complex<utype>> &,         // real,imag in ([0,y],[0,y])
          const utype &) __xprototype;

  template <class utype>
  xtensor6<std::complex<utype>>     // complex xtensor6 of uniforms
  uniform(const xtensor6<std::complex<utype>> &,         // real,imag in ([0,x],[0,y])
          const utype &,
          const utype &) __xprototype;

  template <class utype>
  xtensor6<std::complex<utype>>     // complex xtensor6 of uniforms
  uniform(const xtensor6<std::complex<utype>> &,     // real,imag in ([x1,y1],[x2,y2])
          const utype &,
          const utype &,
          const utype &,
          const utype &) __xprototype;
} // namespace aurostd

// --------------------------------- gaussian function and procedures on xtensor6

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor6<utype>             // gaussian xtensor6<utype> (m,v)
  gaussian(const xtensor6<utype> &a, const utype &m, const utype &v) __xprototype;

  template <class utype>
  xtensor6<utype>           // gaussian xtensor6<utype> (m=0,v)
  gaussian(const xtensor6<utype> &a, const utype &v) __xprototype;

  template <class utype>
  xtensor6<utype>         // gaussian xtensor6<utype> (m=0,v=1)
  gaussian(const xtensor6<utype> &a) __xprototype;

  template <class utype>
  xtensor6<std::complex<utype>>    // complex xtensor6 of gaussians
  gaussian(const xtensor6<std::complex<utype>> &a,    // real,imag in (m1,v1),(m2,v2)
           const utype &m2,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor6<std::complex<utype>>    // complex xtensor6 of gaussians
  gaussian(const xtensor6<std::complex<utype>> &a,      // real,imag in (0,v1),(0,v2)
           const utype &v1,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor6<std::complex<utype>>    // complex xtensor6 of gaussians
  gaussian(const xtensor6<std::complex<utype>> &a,        // real,imag in (0,v),(0,v)
           const utype &v) __xprototype;

  template <class utype>
  xtensor6<std::complex<utype>>    // complex xtensor6 of gaussians
  gaussian(const xtensor6<std::complex<utype>> &a)       // real,imag in (0,1),(0,1)
      __xprototype;
} // namespace aurostd

#endif //__XTENSOR6_CPP

// ---------------------------------- uniform function and procedures on xtensor7

#ifdef __XTENSOR7_CPP

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor7<utype>                 // xtensor7 of uniforms [0,1]
  uniform(const xtensor7<utype> &) __xprototype;

  template <class utype>
  xtensor7<utype>                 // xtensor7 of uniforms [0,y]
  uniform(const xtensor7<utype> &, const utype &) __xprototype;

  template <class utype>
  xtensor7<utype>                 // xtensor7 of uniforms [x,y]
  uniform(const xtensor7<utype> &, const utype &, const utype &) __xprototype;

  template <class utype>
  xtensor7<std::complex<utype>>     // complex xtensor7 of uniforms
  uniform(const xtensor7<std::complex<utype>> &)         // real,imag in ([0,1],[0,1])
      __xprototype;

  template <class utype>
  xtensor7<std::complex<utype>>     // complex xtensor7 of uniforms
  uniform(const xtensor7<std::complex<utype>> &,         // real,imag in ([0,y],[0,y])
          const utype &) __xprototype;

  template <class utype>
  xtensor7<std::complex<utype>>     // complex xtensor7 of uniforms
  uniform(const xtensor7<std::complex<utype>> &,         // real,imag in ([0,x],[0,y])
          const utype &,
          const utype &) __xprototype;

  template <class utype>
  xtensor7<std::complex<utype>>     // complex xtensor7 of uniforms
  uniform(const xtensor7<std::complex<utype>> &,     // real,imag in ([x1,y1],[x2,y2])
          const utype &,
          const utype &,
          const utype &,
          const utype &) __xprototype;
} // namespace aurostd

// --------------------------------- gaussian function and procedures on xtensor7

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor7<utype>             // gaussian xtensor7<utype> (m,v)
  gaussian(const xtensor7<utype> &a, const utype &m, const utype &v) __xprototype;

  template <class utype>
  xtensor7<utype>           // gaussian xtensor7<utype> (m=0,v)
  gaussian(const xtensor7<utype> &a, const utype &v) __xprototype;

  template <class utype>
  xtensor7<utype>         // gaussian xtensor7<utype> (m=0,v=1)
  gaussian(const xtensor7<utype> &a) __xprototype;

  template <class utype>
  xtensor7<std::complex<utype>>    // complex xtensor7 of gaussians
  gaussian(const xtensor7<std::complex<utype>> &a,    // real,imag in (m1,v1),(m2,v2)
           const utype &m2,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor7<std::complex<utype>>    // complex xtensor7 of gaussians
  gaussian(const xtensor7<std::complex<utype>> &a,      // real,imag in (0,v1),(0,v2)
           const utype &v1,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor7<std::complex<utype>>    // complex xtensor7 of gaussians
  gaussian(const xtensor7<std::complex<utype>> &a,        // real,imag in (0,v),(0,v)
           const utype &v) __xprototype;

  template <class utype>
  xtensor7<std::complex<utype>>    // complex xtensor7 of gaussians
  gaussian(const xtensor7<std::complex<utype>> &a)       // real,imag in (0,1),(0,1)
      __xprototype;
} // namespace aurostd

#endif //__XTENSOR7_CPP

// ---------------------------------- uniform function and procedures on xtensor8

#ifdef __XTENSOR8_CPP

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor8<utype>                 // xtensor8 of uniforms [0,1]
  uniform(const xtensor8<utype> &) __xprototype;

  template <class utype>
  xtensor8<utype>                 // xtensor8 of uniforms [0,y]
  uniform(const xtensor8<utype> &, const utype &) __xprototype;

  template <class utype>
  xtensor8<utype>                 // xtensor8 of uniforms [x,y]
  uniform(const xtensor8<utype> &, const utype &, const utype &) __xprototype;

  template <class utype>
  xtensor8<std::complex<utype>>     // complex xtensor8 of uniforms
  uniform(const xtensor8<std::complex<utype>> &)         // real,imag in ([0,1],[0,1])
      __xprototype;

  template <class utype>
  xtensor8<std::complex<utype>>     // complex xtensor8 of uniforms
  uniform(const xtensor8<std::complex<utype>> &,         // real,imag in ([0,y],[0,y])
          const utype &) __xprototype;

  template <class utype>
  xtensor8<std::complex<utype>>     // complex xtensor8 of uniforms
  uniform(const xtensor8<std::complex<utype>> &,         // real,imag in ([0,x],[0,y])
          const utype &,
          const utype &) __xprototype;

  template <class utype>
  xtensor8<std::complex<utype>>     // complex xtensor8 of uniforms
  uniform(const xtensor8<std::complex<utype>> &,     // real,imag in ([x1,y1],[x2,y2])
          const utype &,
          const utype &,
          const utype &,
          const utype &) __xprototype;
} // namespace aurostd

// --------------------------------- gaussian function and procedures on xtensor8

namespace aurostd {
  // namespace aurostd
  template <class utype>
  xtensor8<utype>             // gaussian xtensor8<utype> (m,v)
  gaussian(const xtensor8<utype> &a, const utype &m, const utype &v) __xprototype;

  template <class utype>
  xtensor8<utype>           // gaussian xtensor8<utype> (m=0,v)
  gaussian(const xtensor8<utype> &a, const utype &v) __xprototype;

  template <class utype>
  xtensor8<utype>         // gaussian xtensor8<utype> (m=0,v=1)
  gaussian(const xtensor8<utype> &a) __xprototype;

  template <class utype>
  xtensor8<std::complex<utype>>    // complex xtensor8 of gaussians
  gaussian(const xtensor8<std::complex<utype>> &a,    // real,imag in (m1,v1),(m2,v2)
           const utype &m2,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor8<std::complex<utype>>    // complex xtensor8 of gaussians
  gaussian(const xtensor8<std::complex<utype>> &a,      // real,imag in (0,v1),(0,v2)
           const utype &v1,
           const utype &v2) __xprototype;

  template <class utype>
  xtensor8<std::complex<utype>>    // complex xtensor8 of gaussians
  gaussian(const xtensor8<std::complex<utype>> &a,        // real,imag in (0,v),(0,v)
           const utype &v) __xprototype;

  template <class utype>
  xtensor8<std::complex<utype>>    // complex xtensor8 of gaussians
  gaussian(const xtensor8<std::complex<utype>> &a)       // real,imag in (0,1),(0,1)
      __xprototype;
} // namespace aurostd

#endif //__XTENSOR8_CPP

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XRANDOM_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
