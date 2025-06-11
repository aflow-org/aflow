// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XVECTOR_H_
#define _AUROSTD_XVECTOR_H_
#define _AUROSTD_XVECTOR_DEFAULT_SIZE_ 3


#define BOUNDARY_CONDITIONS_NONE 0
#define BOUNDARY_CONDITIONS_PERIODIC 1    

#define _AUROSTD_XVECTOR_TOLERANCE_IDENTITY_ AUROSTD_IDENTITY_TOL //1.0e-6
#define _AUROSTD_XVECTOR_TOLERANCE_ROUNDOFF_ AUROSTD_ROUNDOFF_TOL //1.0e-6

#define CONV_SHAPE_FULL 0 //CO20190520
#define CONV_SHAPE_SAME 1 //CO20190520
#define CONV_SHAPE_VALID 2 //CO20190520

// -------------------------------------------------------------- class xvector

namespace aurostd {
  template<class utype> class xmatrix;  //forward declaration  //CO20191110
  // namespace aurostd
  template<class utype>
    class xvector {
      public:
        // xvector();                                          // default constructor
        // xvector(int);                                       // default constructor
        xvector(int=3,int=1);                                  // default constructor
        xvector(const std::initializer_list<utype>);           //  initializer_list constructor //HE20220616
        xvector(const xvector<utype>&);                        // copy constructor
        xvector(const xmatrix<utype>&);                        // copy constructor //CO20191110
        // xvector (const xmatrix<utype>&);                    // make a vector of a xmatrix
        xvector<utype>& operator=(const xvector<utype>&);	     // assignment
        xvector<utype>& operator=(const std::initializer_list<utype>); // initializer_list assignment //HE20220616
        // operator xvector<utype>() { return *this;};         // IBM_CPP
        ~xvector();                                            // default destructor
        utype& operator[](int) const;		                       // indicize
        utype& operator()(int) const;		                       // indicize
        utype& operator()(int,bool) const;                     // indicize boundary conditions
        xvector<utype>& operator +=(const xvector<utype>&);
        xvector<utype>& operator +=(const std::initializer_list<utype>); //HE20220616
        xvector<utype>& operator +=(utype); //CO20180409
        xvector<utype>& operator -=(const xvector<utype>&);
        xvector<utype>& operator -=(utype); //CO20180409
        xvector<utype>& operator -=(const std::initializer_list<utype>); //HE20220616
        xvector<utype>& operator *=(utype r); //(const utype& r) //CO20171130 - v*=v[1] doesn't work //CO20180409
        xvector<utype>& operator /=(utype r); //(const utype& r) //CO20171130 - v/=v[1] doesn't work //CO20180409
        // xvector<utype>& operator *=(const xvector<utype>&);
        // xvector<utype>& operator /=(const xvector<utype>&);
        //    friend std::ostream& operator<<<utype>(std::ostream&,const xvector<utype>&);
        //    friend std::ostream& operator< <utype>(std::ostream&,const xvector<utype>&);
        int rows,lrows,urows;  
        bool isfloat,iscomplex;
        // operations
        void set(const utype&);
        void reset(void);
        void clear(void);
        void null(void);  //CO20200731 - to create null vector
        void resize(int=3,int nl=1); //CO20201111
        //HE20220915 START
        // xvector manipulation without copy (pointer manipulation)
        void shift(int new_lrows);
        //HE20220915 END
      private:
        utype *corpus;
        // bool isfloat,iscomplex;
        char size;
        long int vsize;
        typedef typename std::initializer_list<utype>::const_iterator ili; // initializer list iterator //HE20220616

        //NECESSARY PRIVATE CLASS METHODS - START
        void init(); //HE20220515
        void free();  //CO20190808
        void copy(const xvector<utype>& b);  //CO20190808
        void copy(const xmatrix<utype>& b);  //CO20190808
        void copy(const std::initializer_list<utype> l);  //HE20220616
        void refresh(); //CO20190808 - refresh internal properties dependent on lrows, urows, utype

        //NECESSARY END CLASS METHODS - END
    };
}

// ----------------------------------------------------------------------------
// ------------------------- primitives for template<class utype> xvector<utype>

// ----------------------------------------------------------------------------
// ------------------------------------------------- unary and binary operators

namespace aurostd {
  // namespace aurostd
  // template<class utype> std::ostream& operator<<(std::ostream&,const xvector<utype>&);

  template<class utype>
    std::ostream& operator<<(std::ostream&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator+(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator-(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator+(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator-(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator%(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator+(const utype,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator+(const xvector<utype>&,const utype) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator+(const stype,const xvector<utype>&) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator+(const xvector<utype>&,const stype) __xprototype;

  template<class utype> xvector<utype>
    operator-(const utype,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator-(const xvector<utype>&,const utype) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator-(const stype,const xvector<utype>&) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator-(const xvector<utype>&,const stype) __xprototype;

  template<class utype> xvector<utype>
    operator*(const utype,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator*(const xvector<utype>&,const utype) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator*(const stype,const xvector<utype>&) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator*(const xvector<utype>&,const stype) __xprototype;

  template<class utype> xvector<utype>
    operator/(const utype,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator/(const xvector<utype>&,const utype) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator/(const stype,const xvector<utype>&) __xprototype;

  template<class utype,class stype> xvector<utype>
    operator/(const xvector<utype>&,const stype) __xprototype;

  //ME20200329 - real * complex vector
  template<class utype> xvector<xcomplex<utype> >
    operator*(utype, const xvector<xcomplex<utype> >&);

  template<class utype> xvector<xcomplex<utype> >
    operator*(const xvector<xcomplex<utype> >&, utype);

  template<class utype> xvector<xcomplex<utype> >
    operator/(const xvector<xcomplex<utype> >&, utype);

  template<class utype> xvector<utype>
    operator<<(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    operator<<(const xvector<utype>&,const utype) __xprototype;

  template<class utype> xvector<utype>
    operator<<(const utype,const xvector<utype>&) __xprototype;

  template<class utype> utype
    operator*(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> utype
    scalar_product(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    vector_product(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>                
    elementwise_product(const xvector<utype>&,const xvector<utype>&) __xprototype; //SD20220422

  template<class utype> utype
    elements_product(const xvector<utype>&) __xprototype; //SD20220617

  //ME20200327
  template<class utype> xmatrix<utype>
    outer_product(const xvector<utype>&, const xvector<utype>&) __xprototype;

  template<class utype> xvector<char>                    // is xvector > scalar ?
    operator>(const xvector<utype>&,const utype&) __xprototype;

  template<class utype> xvector<char>                    // is xvector < scalar ?
    operator<(const xvector<utype>&,const utype&) __xprototype;

  template<class utype> xvector<char>                    // is xvector == scalar ?
    operator==(const xvector<utype>&,const utype&) __xprototype;

  template<class utype> xvector<char>                    // is xvector > xvector ?
    operator>(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> xvector<char>                    // is xvector < xvector ?
    operator<(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    identical(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;

  template<class utype> bool
    identical(const xvector<xcomplex<utype>>&,const xvector<xcomplex<utype>>&,const utype&) __xprototype;

  template<class utype> bool
    identical(const xvector<utype>&, utype tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype; //DX20210503

  template<class utype> bool
    identical(const xvector<xcomplex<utype>>&, utype tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype; //DX20210503

  template<class utype> bool
    identical(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    identical(const xvector<xcomplex<utype>>&,const xvector<xcomplex<utype>>&) __xprototype;

  template<class utype> bool
    operator==(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    isdifferent(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;

  template<class utype> bool
    isdifferent(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    isequal(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;

  template<class utype> bool
    isequal(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    isequal(const xvector<xcomplex<utype>>&,const xvector<xcomplex<utype>>&) __xprototype;

  template<class utype> bool
    operator!=(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool
    isinteger(const xvector<utype>&,const utype& tol=(utype)0.01) __xprototype; //CO20180409

  template<class utype> bool
    iszero(const xvector<utype>&, double tol=(double)_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype;  //ME20180702 //CO20191201 - 1e-7 seems arbitrary  //CO20230313 - utype clang warnings

}

// ----------------------------------------------------------- xvector example types
namespace aurostd {
  //CO20180419
  template<class utype> xvector<utype> null_xv() __xprototype;  //CO20200731 - friend so it can access refresh()
  template<class utype> xvector<utype> ones_xv(int=3,int=1) __xprototype;
  template<class utype> xvector<utype> box_filter_xv(int window,int lrows=1) __xprototype;
  template<class utype> int gaussian_filter_get_window(utype sigma) __xprototype;  //CO20220627
  template<class utype> xvector<utype> gaussian_filter_xv(utype sigma) __xprototype;  //if you need lrows!=1, use shift()
  template<class utype> xvector<utype> gaussian_filter_xv(utype sigma,int window,int lrows=1) __xprototype; //CO20190419
}

namespace aurostd {

  // CAST ON XVECTORS

  template<class utype> xvector<long double>
    xlongdouble(const xvector<utype>&) __xprototype;

  template<class utype> xvector<double>
    xdouble(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    mod(const xvector<utype>&,utype d) __xprototype;  //CO20200127

  template<class utype> xvector<utype>
    mod_floored(const xvector<utype>&,utype d) __xprototype;  //SD20220117

  template<class utype> xvector<double>
    floor(const xvector<utype>&) __xprototype;

  template<class utype> xvector<double>
    ceil(const xvector<utype>&) __xprototype;

  template<class utype> xvector<double>
    round(const xvector<utype>&) __xprototype;

  template<class utype> xvector<double>
    trunc(const xvector<utype>&) __xprototype;

  template<class utype> xvector<float>
    xfloat(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    pow(const xvector<utype>&,const utype d) __xprototype;  //SD20220324

  template<class utype> xvector<utype>
    pow(const xvector<utype>&,const xvector<utype>&) __xprototype;  //SD20220324

  template<class utype> xvector<long int>
    xlongint(const xvector<utype>&) __xprototype;

  template<class utype> xvector<int>
    xint(const xvector<utype>&) __xprototype;

  template<class utype> xvector<char>
    xchar(const xvector<utype>&) __xprototype;

  template<class utype> vector<utype>
    xvector2vector(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    vector2xvector(const vector<utype>&,int lrows=1) __xprototype; //CO20180409
  template<class utype> xvector<utype>
    vector2xvector(const vector<string>&,int lrows=1) __xprototype; //CO20180409

  //[SD20220512 - OBSOLETE]xvector<double> xvectorint2double(const xvector<int>&); //CO20180515
  //[SD20220512 - OBSOLETE]xvector<int> xvectordouble2int(const xvector<double>&,bool check_int=true); //CO20180515

  template<class utype1,class utype2> xvector<utype2>
    xvector2utype(const xvector<utype1>& a) __xprototype; //SD20220512

  // OPERATIONS ON XVECTORS

  template<class utype> void
    set(xvector<utype>&, const utype&) __xprototype;

  template<class utype> void
    reset(xvector<utype>&) __xprototype;

  template<class utype> void
    clear(xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    vabs(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    abs(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    sign(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    nint(const xvector<utype>&) __xprototype;

  template<class utype> utype
    modulus(const xvector<utype>&) __xprototype;

  template<class utype> utype
    modulussquare(const xvector<utype>&) __xprototype;

  template<class utype> utype
    modulus2(const xvector<utype>&) __xprototype;

  template<class utype> utype
    sum(const xvector<utype>&) __xprototype;

  template<class utype> utype
    min(const xvector<utype>&) __xprototype;

  template<class utype> utype
    min(const xvector<utype>&,int&) __xprototype;

  template<class utype> int
    mini(const xvector<utype>&) __xprototype;

  template<class utype> utype
    max(const xvector<utype>&) __xprototype;

  template<class utype> utype
    max(const xvector<utype>&,int&) __xprototype;

  template<class utype> int
    maxi(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>                  // clear values too small
    roundoff(const xvector<utype>&,utype tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_ROUNDOFF_) __xprototype; // claar values too small //CO20180409

  void GCD(const xvector<int>&,int&);                         // get GCD of xvector //CO20180409 //CO20191201
  void GCD(const xvector<int>& va,const xvector<int>& vb,xvector<int>& vgcd); //CO20191201
  void GCD(const xvector<int>& va,const xvector<int>& vb,xvector<int>& vgcd,xvector<int>& vx,xvector<int>& vy); //CO20191201
  int LCM(const xvector<int>&);                         // get LCM of xvector //CO20180520

  //DX20191125 [OBSOLETE] template<class utype> xvector<utype>                                       // simply divide by GCD (useful for compounds) //CO20180409
  //DX20191125 [OBSOLETE]   reduceByGCD(const xvector<utype>& in_V,                                  // simply divide by GCD (useful for compounds) //CO20180409
  //DX20191125 [OBSOLETE]       const utype& tol=_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype; // simply divide by GCD (useful for compounds) //CO20180409

  template<class utype> 
    void reduceByGCD(const xvector<utype>& in_V,                    // simply divide by GCD (useful for compounds) //CO20180409
        xvector<utype>& out_V,                                        // simply divide by GCD (useful for compounds) //CO20180409
        utype tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype; // simply divide by GCD (useful for compounds) //CO20180409 //CO20191201

  void GCD(const vector<int>&,int&);                          // get GCD of vector //DX20191125 (modeled after CO's xvector version) //CO20191201
  int LCM(const vector<int>&);                          // get LCM of vector //DX20191125 (modeled after CO's xvector version)

  template<class utype> 
    void reduceByGCD(const vector<utype>& in_V,                     // simply divide by GCD (useful for compounds) //DX20191125 (modeled after CO's xvector version)
        vector<utype>& out_V,                                         // simply divide by GCD (useful for compounds) //DX20191125 (modeled after CO's xvector version)
        utype tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype; // simply divide by GCD (useful for compounds) //DX20191125 (modeled after CO's xvector version)  //CO20191201

  void GCD(const deque<int>&,int& gcd);                           // get GCD of deque //DX20191125 (modeled after CO's xvector version)  //CO20191201
  int LCM(const deque<int>&);                           // get LCM of deque //DX20191125 (modeled after CO's xvector version)

  template<class utype> 
    void reduceByGCD(const deque<utype>& in_V,                      // simply divide by GCD (useful for compounds) //DX20191125 (modeled after CO's xvector version)
        deque<utype>& out_V,                                          // simply divide by GCD (useful for compounds) //DX20191125 (modeled after CO's xvector version)
        utype tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype; // simply divide by GCD (useful for compounds) //DX20191125 (modeled after CO's xvector version)  //CO20191201

  template<class utype> xvector<utype> 
    normalizeSumToOne(const xvector<utype>& in_V,
        const utype& tol=(utype)_AUROSTD_XVECTOR_TOLERANCE_ROUNDOFF_) __xprototype; //CO20180801

  template<class utype> void                                   // swap
    swap(xvector<utype>&,const int&,const int&) __xprototype;  // swap

  // Complex operations
  template<class utype> xvector<utype>
    conj(const xvector<utype>&) __xprototype;

  // EXPONENTIAL OPERATIONS

  template<class utype> xvector<utype>
    exp(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    log(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    log10(const xvector<utype>&) __xprototype;

  // TRIDIMENSIONAL OPERATIONS

  template<class utype> utype
    distance(const xvector<utype>&,const xvector<utype>&) __xprototype;

  // TRIGONOMETRIC OPERATIONS

  template<class utype> xvector<utype>
    sin(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cos(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    asin(const xvector<utype>&) __xprototype;
  

  template<class utype> xvector<utype>
    acos(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    tan(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    atan(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    sec(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cosec(const xvector<utype>&) __xprototype;

  // HYPERBOLIC TRIGONOMETRIC OPERATIONS

  template<class utype> xvector<utype>
    sinh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    asinh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cosh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    acosh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    tanh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    atanh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cotanh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    acotanh(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    sech(const xvector<utype>&) __xprototype;

  template<class utype> xvector<utype>
    cosech(const xvector<utype>&) __xprototype;

  // TRIGONOMETRIC OPERATIONS BETWEEN TWO VECTORS

  template<class utype> double   // cos of angle between two vectors
    cos(const xvector<utype>&,const xvector<utype>&) __xprototype;
  template<class utype> double   // cos of angle between two vectors
    getcos(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> double   // sin of angle between two vectors
    sin(const xvector<utype>&,const xvector<utype>&) __xprototype;
  template<class utype> double   // sin of angle between two vectors
    getsin(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> double   // angle between two vectors in radiants
    angle(const xvector<utype>&,const xvector<utype>&) __xprototype;
  template<class utype> double   // angle between two vectors in radiants
    getangle(const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> double   // angle between three vectors in radiants
    angle(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;
  template<class utype> double   // angle between three  vectors in radiants
    getangle(const xvector<utype>&,const xvector<utype>&,const xvector<utype>&) __xprototype;

  template<class utype> bool //determine if two vectors are collinear //CO20180409
    isCollinear(const xvector<utype>& v0,const xvector<utype>& v1,const utype& tol=_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_) __xprototype;

  template<class utype> double // area encapsulated by ordered points (on the same plane) //HE20210511
    areaPointsOnPlane(const vector<xvector<utype> >& points) __xprototype;
  template<class utype> double // volume defined by points, facets and all outward or inward facing normals (needs to be consistent) //HE20210511
    volume(const vector<xvector<utype> >& points, const vector<vector<uint> >& facets, const vector<xvector<utype> >& normals) __xprototype;
  double // volume defined by points and ordered facets //HE20210511
    volume(const vector<xvector<int> >& points, const vector<vector<uint> >& facets, bool convex = false);
  template<class utype> double // volume defined by points and ordered facets //HE20210511
    volume(const vector<xvector<utype> >& points, const vector<vector<uint> >& facets, bool convex = false) __xprototype;

  template<class utype> xvector<utype> //get centroid of data points //CO20180409
    getCentroid(const vector<xvector<utype> >& points) __xprototype;
  template<class utype> xvector<utype> //get centroid of data points with weights //DX20200728
    getCentroid(const vector<xvector<utype> >& points, const vector<utype>& weights) __xprototype;

  template<class utype> xvector<double> //get centroid of data points //CO20180409
    getCentroidPBC(const vector<xvector<utype> >& points, const xmatrix<utype>& lattice) __xprototype;
  template<class utype> xvector<double> //get centroid of data points with weights //DX20200728
    getCentroidPBC(const vector<xvector<utype> >& points, const vector<utype>& weights, const xmatrix<utype>& lattice) __xprototype;

  // TRIGONOMETRIC OPERATIONS BETWEEN TWO GENERAL VECTORS //CO20180409
  template<class utype> xvector<double> 
    getGeneralAngles(const xvector<utype>& vec,const utype& tol=_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_); //CO20180409

  template<class utype> double
    getGeneralAngle(const xvector<utype>& _vec,int _i,const utype& tol=_AUROSTD_XVECTOR_TOLERANCE_IDENTITY_); //CO20180409

  template<class utype> xvector<utype>
    getGeneralNormal(const vector<xvector<utype> >& _directive_vectors); //CO20180409

  template<class utype> xvector<utype>
    pointLineIntersection(const xvector<utype>& a,const xvector<utype>& n,const xvector<utype>& p); //CO20180520

  template<class utype> bool
    linePlaneIntersect(const xvector<utype>& p0,const xvector<utype>& n,const xvector<utype>& l0, const xvector<utype>& l,double& d,xvector<utype>& intersection); //CO20180520

  template<class utype> xvector<utype> getVectorProjection(const xvector<utype>& b, const xvector<utype>& a);  //ME20200511
  template<class utype> xvector<utype> getModeratedVectorProjection(const xvector<utype> c, const xvector<utype>& b, const xvector<utype>& a);  //ME20200511

  // grid generation
  xvector<double> linspace(const double start, const double stop, const int n);
  xvector<double> linspace(const double start, const double stop, const double n);

  // SIMPLE SORT ROUTINES
  template<class utype> void // WRAP TO SHELL SHORT
    sort(xvector<utype>& a) __xprototype;

  template<class utype> xvector<utype> // SHELLSORT
    shellsort(const xvector<utype>& a) __xprototype;

  template<class utype> xvector<utype> // HEAPSORT
    heapsort(const xvector<utype>& a) __xprototype;

  template<class utype> xvector<utype>  // QUICKSORT
    quicksort(const xvector<utype>&) __xprototype;

  template<class utype1, class utype2> void  // QUICKSORT
    quicksort2(unsigned long n, xvector<utype1>&, xvector<utype2>&) __xprototype;

  template<class utype1, class utype2, class utype3> void   // QUICKSORT
    quicksort3(unsigned long n, xvector<utype1>&, xvector<utype2>&, xvector<utype3>&)
    __xprototype;

  template<class utype1, class utype2, class utype3> void   // QUICKSORT
    quicksort4(unsigned long n, xvector<utype1>&, xvector<utype2>&, xvector<utype3>&)
    __xprototype;

  template<class utype> xvector<utype>             // SSORT shortcup for SHELLSORT
    ssort(const xvector<utype>& a) { return shellsort(a);}

  template<class utype> xvector<utype>             // HPSORT shortcup for HEAPSORT  
    hpsort(const xvector<utype>& a) { return heapsort(a);}

  //  template<class utype> xvector<utype>             // QSORT shortcup for QUICKSORT  
  //  sort(const xvector<utype>& a) { return heapsort(a);}

  template<class utype> xvector<utype>                 // SORT shortcup for HPSORT  
    qsort(const xvector<utype>& a) { return quicksort(a);}

  template<class utype1, class utype2> void  // QUICKSORT
    sort2(unsigned long n, xvector<utype1>&a, xvector<utype2>&b) {
      quicksort2(n,a,b);}

  template<class utype1, class utype2, class utype3> void  // QUICKSORT
    sort3(unsigned long n, xvector<utype1>&a, xvector<utype2>&b, xvector<utype3>&c) {
      quicksort3(n,a,b,c);}

  template<class utype1, class utype2, class utype3, class utype4> void  // QUICKSORT
    sort4(unsigned long n, xvector<utype1>&a, xvector<utype2>&b, xvector<utype3>&c, xvector<utype4>&d) {
      quicksort4(n,a,b,c,d);}
}

namespace aurostd { //CO20190419
  template<class utype>
    class compareVecElement {
      public:
        compareVecElement(uint ind=0,bool ascending=true);
        bool operator() (const vector<utype>& a,const vector<utype>& b);
        bool operator() (const xvector<utype>& a,const xvector<utype>& b);
      private:
        uint m_uindex_sort; //keep in memory so we don't rely on many conversions per sort
        int m_iindex_sort;  //keep in memory so we don't rely on many conversions per sort
        bool m_ascending_sort;  //m_ascending_sort==true is ascending sort, ==false is descending sort, NB this can be different than using rbegin()/rend(), see sort(ids...) in XRD analysis in aflow_pflow_funcs.cpp 
    };
  template<class utype> bool compareVecElements(const vector<utype>& a,const vector<utype>& b);
  template<class utype> bool compareXVecElements(const aurostd::xvector<utype>& a,const aurostd::xvector<utype>& b);
}

namespace aurostd { //CO20190419
  //SOME STATS STUFF
  template<class utype> utype mean(const xvector<utype>& a); //CO20190520
  template<class utype> utype meanWeighted(const xvector<utype>& a,const xvector<utype>& weights); //CO20190520
  template<class utype> utype meanWeighted(const xvector<utype>& a,const xvector<utype>& weights,utype& sum_weights); //CO20190520
  template<class utype> utype var(const xvector<utype>& a,int ddof=1); //CO20190520
  template<class utype> utype stddev(const xvector<utype>& a,int ddof=1); //CO20190520
  template<class utype> utype mode(const xvector<utype>& a); //CO20190520
  template<class utype> utype correlation_Pearson_fast(const xvector<utype>& a,const xvector<utype>& b,int ddof=1); //CO20190520
  template<class utype> utype correlation_Pearson_fast(const xvector<utype>& a,utype mean_a,utype stddev_a, const xvector<utype>& b,utype mean_b,utype stddev_b,int ddof=1); //CO20190520
  template<class utype> utype correlation_Pearson_slow(const xvector<utype>& a,const xvector<utype>& b); //CO20190520
  template<class utype> void getQuartiles(const xvector<utype>& _a,utype& q1,utype& q2,utype& q3);  //CO20171202
  template<class utype> utype getMAD(const xvector<utype>& _a,utype median=(utype)AUROSTD_NAN);   //CO20171202, absolute deviation around the median (MAD)
  template<class utype> xvector<utype> convolution(const xvector<utype>& signal_input,const xvector<utype>& response_input,int SHAPE=CONV_SHAPE_FULL); //CO20190419
  template<class utype> xvector<utype> convolution(const xvector<utype>& signal_input,const xvector<utype>& response_input,vector<uint>& sum_counts,int SHAPE=CONV_SHAPE_FULL); //CO20190419
  template<class utype> xvector<utype> moving_average(const xvector<utype>& signal_input,int window); //CO20190419
}

namespace aurostd { //CO20190620
  //signal processing
  template<class utype> vector<int> getPeaks(const xvector<utype>& signal_input,uint smoothing_iterations=4,uint avg_window=4,int width_maximum=1,double significance_multiplier=1.0);  //CO20190620
  template<class utype> vector<int> getPeaks(const xvector<utype>& signal_input,xvector<utype>& signal_smooth,uint smoothing_iterations=4,uint avg_window=4,int width_maximum=1,double significance_multiplier=1.0);  //CO20190620

  template<class utype> vector<xvector<double>> histogram(const xvector<utype>& data, uint bins, double minimum_data, double maximum_data);
  template<class utype> vector<xvector<double>> histogram(const xvector<utype>& data, uint bins); //AZ20230204
  template<class utype> vector<xvector<double>> histogram(const xvector<utype>& data, uint max_bins, uint bin_mode); //AZ20230204
														     
}

namespace aurostd {
  // differenitation
  xvector<double> diffSG(const xvector<double> &f, double dx);//AS20210901
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XVECTOR_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

