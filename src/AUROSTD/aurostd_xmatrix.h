// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XMATRIX_H_
#define _AUROSTD_XMATRIX_H_

#define _AUROSTD_XMATRIX_DEFAULT_SIZE_ 3
//#define _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_ (1.0e-6)
//#define _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_ (1.0e-6)

//double _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_ = 1.0e-6;
//double _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_=1.0e-6;

#ifndef _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_
#define _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_ AUROSTD_ROUNDOFF_TOL //1.0e-6 //DX20171025
#endif

#ifndef _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_
#define _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_ AUROSTD_IDENTITY_TOL //1.0e-6 //DX20171025
#endif

const int _MAX_ITS=100; // maximum iterations in SVDcmp()

// ----------------------------------------------------------------------------
// -------------------------------------------------------------- class xmatrix

namespace aurostd {
  // namespace aurostd
  template<class utype>
    class xmatrix  {
      public:
        // xmatrix ();                                                      // default constructor
        // xmatrix (int=3);                                                 // default constructor
        xmatrix (int=3,int=3,int=1,int=1);                                  // default constructor
        // xmatrix (int,int,int,int);                                       // default constructor
        xmatrix(const std::initializer_list<std::initializer_list<utype>>); // initializer_list constructor //HE20220616
        xmatrix (int,int,utype*);                                           // copy constructor
        xmatrix (const xmatrix<utype>&);                                    // copy constructor
        xmatrix (const xvector<utype>&);                                    // copy constructor
        xmatrix<utype>& operator=(const xmatrix<utype>&);                   // assignment
        xmatrix<utype>& operator=(const std::initializer_list<std::initializer_list<utype>>); // initializer_list assignment //HE20220616
        ~xmatrix ();                                                        // default destructor
        utype* operator[](int) const;                                       // indicize i,j
        xvector<utype> operator()(int) const;                               // indicize i
        xvector<utype> getcol(int) const;                                   // return column i  //CO20191110
        xvector<utype> getdiag(int k=0,int _lrows=1) const;                 // return diagonal, k=0 is main diagonal, k>0 is above diagonal //CO20191201
        xvector<utype> getxvec() const;                                      // updated function to return a vector from a 1-d matrix //AZ20220704
        xvector<utype> getxvec(int lrow,int urow,int lcol,int ucol, int lrows_out=1) const;   // return a vector from a matrix //AZ20220704
        void getxmatInPlace(xmatrix<utype>& mat_out,int lrow,int urow,int lcol,int ucol,int lrows_out=1,int lcols_out=1) const; // return submatrix as xmatrix //CO20191110
        void getxvecInPlace(xvector<utype>& xv_out,int lrow,int urow,int lcol,int ucol, int lrows_out=1) const; // return xvectorInPlace from xmatrix 
        xmatrix<utype> getxmat(int lrow,int urow,int lcol,int ucol,int lrows_out=1,int lcols_out=1) const; // return submatrix as xmatrix //CO20191110
        void setrow(const xvector<utype>& row,int irow=1);                  // set row of matrix //CO20190808
        void setcol(const xvector<utype>& col,int icol=1);                  // set col of matrix //CO20190808
        void setmat(const xmatrix<utype>& mat,int irow=1,int icol=1);       // set submat  //CO20190808
        void setmat(const xvector<utype>& xv,int icol=1,bool col=true);     // set submat  //CO20190808
        utype& operator()(int,int) const;                                   // indicize i,j
        utype& operator()(int,int,bool) const;                              // indicize boundary conditions
        // xmatrix operator()(int,int,int,int) const;                       // indicize i1,j1,i2,j2
        // math operators
        xmatrix<utype>& operator +=(const xmatrix<utype>&);
        xmatrix<utype>& operator +=(const std::initializer_list<std::initializer_list<utype>>); //HE20220616
        xmatrix<utype>& operator -=(const xmatrix<utype>&);
        xmatrix<utype>& operator -=(const std::initializer_list<std::initializer_list<utype>>); //HE20220616
        xmatrix<utype>& operator *=(const xmatrix<utype>&);
        xmatrix<utype>& operator *=(utype r); //CO20190911
        // xmatrix<utype>& operator /=(const xmatrix<utype>&);
        xmatrix<utype>& operator /=(utype r); //CO20190911
        xmatrix<utype>& operator /=(const xmatrix<utype>&); //CO20191201 - right matrix division
        // std::ostream operator
        //      friend std::ostream& operator<<<utype>(std::ostream&,const xmatrix<utype>&);
        // friend std::ostream& operator<<<utype>(std::ostream&,const xmatrix<utype>&);
        // friend void balanc(const xmatrix<utype>&); // balanced with eigenvalues
        int rows,lrows,urows;                                  // rows dimensions
        int cols,lcols,ucols;                                  // cols dimensions
        bool issquare,isfloat,iscomplex;
        // operations
        void set(const utype&);
        void reset(void);
        void clear(void);
        //HE20220915 START
        // xmatrix manipulation without copy (pointer manipulation)
        void shift(int new_lrows, int new_lcols);
        void shiftrow(int new_lrows);
        void shiftcol(int new_lcols);
        void reshape(uint new_rows, uint new_cols=0);
        void reshape(int new_urows, int new_ucols, int new_lrows, int new_lcols);
        //HE20220915 END
      private:
        //  friend bool operator==(const xmatrix<utype>&,const xmatrix<utype>&);
        utype** corpus;
        //  bool isfloat,iscomplex;
        char size;
        long int msize;
        typedef typename std::initializer_list<utype>::const_iterator ili; // initializer list iterator //HE20220616
        typedef typename std::initializer_list<std::initializer_list<utype>>::const_iterator il2i; // nested initializer list iterator //HE20220616
        //NECESSARY PRIVATE CLASS METHODS - START
        void init();  //HE20220615
        void free();  //CO20190808
        void copy(const xmatrix<utype>& b);  //CO20190808
        void copy(const xvector<utype>& b);  //CO20190808
        void copy(const std::initializer_list<std::initializer_list<utype>> ll); //HE20220616
        void refresh(); //CO20190808 - refresh internal properties dependent on lrows, urows, utype
        //NECESSARY END CLASS METHODS - END
    };
}

// ----------------------------------------------------------------------------
// ------------------------- primitives for template<class utype> xmatrix<utype>

// ------------------------------------------------- unary and binary operators
namespace aurostd {
  template<class utype> xmatrix<utype>
    operator+(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    operator-(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>                                   // matrix*matrix
    operator*(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype; // matrix*matrix

  //ME20190814
  template<class utype> xmatrix<xcomplex<utype> >                        // matrix*matrix
    operator*(const xmatrix<utype>&,const xmatrix<xcomplex<utype> >&) __xprototype; // matrix*matrix

  template<class utype> xvector<utype>                                   // matrix*vector
    operator*(const xmatrix<utype>&,const xvector<utype>&) __xprototype; // matrix*vector

  template<class utype> xvector<utype>                                   // vector_row*matrix
    operator*(const xvector<utype>&,const xmatrix<utype>&) __xprototype; // vector_row*matrix

  //ME20200330 - Multiplication of a real matrix with a complex vector
  template<class utype> xvector<xcomplex<utype> >
    operator*(const xmatrix<utype>&, const xvector<xcomplex<utype> >&) __xprototype;

  template<class utype> xmatrix<utype>                                   // scalar+matrix
    operator+(const utype,const xmatrix<utype>&) __xprototype;           // scalar+matrix

  template<class utype> xmatrix<utype>                                   // matrix+scalar
    operator+(const xmatrix<utype>&,const utype) __xprototype;           // matrix+scalar

  template<class utype> xmatrix<utype>                                   // scalar*matrix
    operator*(const utype,const xmatrix<utype>&) __xprototype;           // scalar*matrix

  template<class utype> xmatrix<utype>                                   // matrix*scalar
    operator*(const xmatrix<utype>&,const utype) __xprototype;           // matrix*scalar

  //ME20200329 - real * complex matrix
  template<class utype> xmatrix<xcomplex<utype> >
    operator*(const xmatrix<xcomplex<utype> >&, utype) __xprototype;

  template<class utype> xmatrix<xcomplex<utype> >
    operator*(utype, const xmatrix<xcomplex<utype> >&) __xprototype;

  template<class utype> xmatrix<xcomplex<utype> >
    operator/(const xmatrix<xcomplex<utype> >&, utype) __xprototype;

  template<class utype> xmatrix<utype>
    operator/(const xmatrix<utype>&,const utype) __xprototype;

  template<class utype> xmatrix<utype>
    operator/(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;  //CO20191201 - right matrix division

  template<class utype>
    std::ostream& operator<<(std::ostream&,const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    operator+(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    operator-(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&,const char& _mode_) __xprototype;

  template<class utype> bool
    identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;

  template<class utype> bool
    rel_identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;

  template<class utype> bool
    abs_identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;

  template<class utype> bool
    identical(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    operator==(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    isdifferent(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;

  template<class utype> bool
    isdifferent(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    isequal(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;

  template<class utype> bool
    isequal(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    operator!=(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;

  //ME20190814 - xcomplex versions
  template<class utype> bool
    identical(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&,const utype&,const char& _mode_) __xprototype;

  template<class utype> bool
    identical(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&,const utype&) __xprototype;

  template<class utype> bool
    rel_identical(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&,const utype&) __xprototype;

  template<class utype> bool
    abs_identical(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&,const utype&) __xprototype;

  template<class utype> bool
    identical(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&) __xprototype;

  template<class utype> bool
    operator==(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&) __xprototype;

  template<class utype> bool
    isdifferent(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&,const utype&) __xprototype;

  template<class utype> bool
    isdifferent(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&) __xprototype;

  template<class utype> bool
    isequal(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&,const utype&) __xprototype;

  template<class utype> bool
    isequal(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&) __xprototype;

  template<class utype> bool
    operator!=(const xmatrix<xcomplex<utype> >&,const xmatrix<xcomplex<utype> >&) __xprototype;

  //CO START
  template<class utype> bool
    isidentity(const xmatrix<utype>&) __xprototype;
  //CO END

  template<class utype> bool
    isdiagonal(const xmatrix<utype>&,const utype& _eps_=(utype)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) __xprototype; //DX20171025

  template<class utype> bool
    isinteger(const xmatrix<utype>&,const utype& tol=(utype)0.01) __xprototype;

  template<class utype> bool
    issymmetric(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    isantisymmetric(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    ishermitian(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    isantihermitian(const xmatrix<utype>&) __xprototype;

}

// ----------------------------------------------------------- xmatrix constrtuction
namespace aurostd {
  template<class utype> xmatrix<utype>
    reshape(const xvector<utype>& vec, int rows, int cols=0) __xprototype;  // ME20210505

  template<class utype> xmatrix<utype>
    reshape(const xvector<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    reshape_cols(const xvector<utype>& vec) __xprototype;

  template<class utype> xmatrix<utype>
    reshape_rows(const xvector<utype>& vec) __xprototype;

  template<class utype> xmatrix<utype>
    hstack(const vector<xvector<utype>>& vec_list) __xprototype;

  template<class utype> xmatrix<utype>
    vstack(const vector<xvector<utype>>& vec_list) __xprototype;

}

// ----------------------------------------------------------- xmatrix example types
namespace aurostd {
  //20171008 - CO
  template<class utype> xmatrix<utype>
    eye(int=3,int=AUROSTD_MAX_INT,int=1,int=1) __xprototype; //CO20190520
  template<class utype> xmatrix<utype>
    ones_xm(int=3,int=3,int=1,int=1) __xprototype; //CO20190520
}

// ----------------------------------------------------------- xmatrix cast
namespace aurostd {
  template<class utype> xmatrix<long double>
    xlongdouble(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<double>
    xdouble(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<float>
    xfloat(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<long int>
    xlongint(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<int>
    xint(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<char>
    xchar(const xmatrix<utype>&) __xprototype;

  template<class utype> vector<vector<utype> >
    xmatrix2vectorvector(const xmatrix<utype>& xmat) __xprototype;

  template<class utype> xmatrix<utype>
    vectorvector2xmatrix(const vector<vector<utype> >&, int lrows=1, int lcols=1) __xprototype;
  template<class utype> xmatrix<utype>
    vectorvector2xmatrix(const vector<vector<string> >&, int lrows=1, int lcols=1) __xprototype;

  template<class utype> xmatrix<double> 
    xmatrixutype2double(const xmatrix<utype>& a); //CO20191201
  template<class utype> xmatrix<utype> 
    xmatrixdouble2utype(const xmatrix<double>& a,bool check_int=true); //CO20191201
}

// ----------------------------------------------------------- xmatrix functions
namespace aurostd {
  template<class utype> void
    reset(xmatrix<utype>&) __xprototype;

  template<class utype> void
    clear(xmatrix<utype>&) __xprototype;

  template<class utype> void
    set(xmatrix<utype>&,const utype&) __xprototype;

  //  template<class utype>                        
  //  xvector<utype> vector(const xmatrix<utype>&) __xprototype;

  template<class utype> utype
    det(const xmatrix<utype>&) __xprototype;

  template<class utype> utype
    determinant(const xmatrix<utype>&) __xprototype;

  template<class utype> void
    submatrixInPlace(const xmatrix<utype>&,xmatrix<utype>&,int,int) __xprototype; //CO20191201 - IN PLACE

  template<class utype> xmatrix<utype>
    submatrix(const xmatrix<utype>&,int,int) __xprototype;

  template<class utype> utype
    minordet(const xmatrix<utype>&,const int&,const int&) __xprototype;

  template<class utype> utype
    minordeterminant(const xmatrix<utype>&,const int&,const int&) __xprototype;

  template<class utype> void
    adjointInPlace(const xmatrix<utype>&,xmatrix<utype>&) __xprototype;  //CO20191201

  template<class utype> xmatrix<utype>
    adjoint(const xmatrix<utype>&) __xprototype;  //CO20191201

  template<class utype> xmatrix<utype>
    inverseByAdjoint(const xmatrix<utype>&) __xprototype;  //CO20191201

  template<class utype> xmatrix<utype>
    inverse(const xmatrix<utype>&) __xprototype;

  template<class utype> bool
    isNonInvertible(const xmatrix<utype>&, xmatrix<utype>&) __xprototype; // RETURN ERROR if non invertible

  template<class utype> xmatrix<utype>
    inverseByQR(const xmatrix<utype>&) __xprototype;  //SD20220426

  template<class utype> void
    LUPDecomposition(const xmatrix<utype>& A, xmatrix<double>& LU, xmatrix<double>& P, utype tol=(utype)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) __xprototype; //SD20220426  //CO20230313 - utype clang warnings
  template<class utype> void
    LUPDecomposition(const xmatrix<utype>& A, xmatrix<double>& L, xmatrix<double>& U, xmatrix<double>& P, utype tol=(utype)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) __xprototype; //SD20220426 //CO20230313 - utype clang warnings

  template<class utype> xmatrix<utype>
    inverseByLUP(const xmatrix<utype>&) __xprototype;  //SD20220426

  template<class utype> utype
    condition_number(const xmatrix<utype>&) __xprototype; //SD20220425

  template<class utype> xmatrix<utype>                  // clear values too small
    roundoff(const xmatrix<utype>&,utype) __xprototype; // claar values too small

  template<class utype> xmatrix<utype>                  // clear values too small
    roundoff(const xmatrix<utype>&) __xprototype;       // clear values too small

  template<class utype> utype
    modulus(const xmatrix<utype>&) __xprototype;  //CO20191110

  template<class utype> utype
    modulussquare(const xmatrix<utype>&) __xprototype;  //CO20191110

  template<class utype> utype
    modulus2(const xmatrix<utype>&) __xprototype; //CO20191110

  template<class utype> utype
    sum(const xmatrix<utype>&) __xprototype;

  template<class utype> xvector<utype>
    sum_column(const xmatrix<utype>&) __xprototype;

  template<class utype> xvector<utype>
    mean_column(const xmatrix<utype>&) __xprototype;

  template<class utype>
    xvector<utype> sum_row(const xmatrix<utype>&) __xprototype;

  template<class utype>
    xvector<utype> mean_row(const xmatrix<utype>&) __xprototype;

  template<class utype> utype
    min(const xmatrix<utype>&) __xprototype;

  template<class utype> utype
    min(const xmatrix<utype>&,int&,int&) __xprototype;

  template<class utype> utype
    max(const xmatrix<utype>&) __xprototype;

  template<class utype> utype
    max(const xmatrix<utype>&,int&,int&) __xprototype;

  //DX20180115 [OBSOLETE] template<class utype> double
  //DX20180115 [OBSOLETE]  trace(const xmatrix<utype>&) __xprototype;
  template<class utype> utype                     //DX20180115 - double to utype
    trace(const xmatrix<utype>&) __xprototype;    //DX20180115 - double to utype

  template<class utype> xmatrix<utype>
    HadamardProduct(const xmatrix<utype>&, const xmatrix<utype>&) __xprototype;  //SD20220422

  template<class utype> xmatrix<utype>
    KroneckerProduct(const xmatrix<utype>&, const xmatrix<utype>&) __xprototype;  //ME20180614

  template<class utype> xmatrix<utype>
    shift_left(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    shift_right(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    shift_up(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    shift_down(const xmatrix<utype>&) __xprototype;

  template<class utype> void                                        // swap_cols
    swap_cols(xmatrix<utype>&,const int&,const int&) __xprototype;  // swap_cols

  template<class utype> void                                           // swap_columns
    swap_columns(xmatrix<utype>&,const int&,const int&) __xprototype;  // swap_columns

  template<class utype> void                                           // swap_rows
    swap_rows(xmatrix<utype>&,const int&,const int&) __xprototype;    // swap_rows

    template<class utype> xmatrix<utype>
    identity(const xmatrix<utype>&) __xprototype;

  // template<class utype> xmatrix<utype>
  //  identity(const int&,const int&,utype) __xprototype;

  template<class utype> xmatrix<utype>
    identity(const utype&, int, int=1);  //ME20200123
  //identity(const utype&,const int&,const int&) __xprototype; // OBSOLETE ME20200123

  template<class utype> xmatrix<utype>  // conjugate
    conj(const xmatrix<utype>&) __xprototype; //ME20180904

  template<class utype> void  // transpose
    traspSquareInPlace(xmatrix<utype>&,xmatrix<utype>&, bool=true) __xprototype; //CO20191201

  template<class utype> void  // transpose
    traspInPlace(const xmatrix<utype>&,xmatrix<utype>&, bool=true) __xprototype; //CO20191201

  template<class utype> void  // transpose
    traspInPlace(xmatrix<utype>&, bool=true) __xprototype; //CO20191201

  template<class utype> xmatrix<utype>  // transpose
    trasp(const xmatrix<utype>&, bool=true) __xprototype; // 25 january 2000, ME20190813 - modified, bool is conjugate (true == trasp with conjugate)

  template<class utype> xmatrix<utype>  // transpose
    trasp(const xvector<utype>&, bool=true) __xprototype; // 5 febrary 2000, ME20190813 - modified, bool is conjugate (true == trasp with conjugate)

  // Mathematical operations

  template<class utype> xmatrix<utype>       // 10 Aug 2007
    abs(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // 10 Oct 2007
    mabs(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Sept 2008
    sign(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Sept 2008
    nint(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Oct 2014
    floor(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Oct 2014
    trunc(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Oct 2014
    round(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>       // Oct 2014
    ceil(const xmatrix<utype>&) __xprototype;

  template<class utype> xmatrix<utype>
    exp(const xmatrix<utype>&) __xprototype; // 31 july 1999, fix 10 Aug 2007

  template<class utype> xmatrix<utype>
    log(const xmatrix<utype>&) __xprototype;// MISSING, probably never

  template<class utype> xmatrix<utype>
    log10(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    sin(const xmatrix<utype>&) __xprototype; // 10 august 1999

  template<class utype> xmatrix<utype>
    asin(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    cos(const xmatrix<utype>&) __xprototype; // 10 august 1999

  template<class utype> xmatrix<utype>
    sec(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    cosec(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    asec(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    acosec(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    acos(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    tan(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    atan(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    sinh(const xmatrix<utype>&) __xprototype;  // 12 august 1999

  template<class utype> xmatrix<utype>
    asinh(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    cosh(const xmatrix<utype>&) __xprototype;   // 12 august 1999

  template<class utype> xmatrix<utype>
    acosh(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    tanh(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype> xmatrix<utype>
    atanh(const xmatrix<utype>&) __xprototype;   // MISSING

  template<class utype>
    xvector<utype> LinearLeastSquares(const xmatrix<utype>&, const xvector<utype>&) __xprototype; //SD20220426

  template<class utype>                          // 25 april 2007
    void GaussJordan(xmatrix<utype>&, xmatrix<utype>&) __xprototype;

  // least square stuff aurostd adaptation of nrecipes
  template<class utype>                          // 1 August 2014
    void gaussj(xmatrix<utype>& a, int n, xmatrix<utype>& b, int m);

  template<class utype>                          // 1 August 2014
    void covsrt(xmatrix<utype>& covar, xvector<int> ia, int mfit);

  template<class utype>                          // 1 August 2014
    void lfit(xvector<utype> x, xvector<utype> y, xvector<utype> sig, 
        xvector<utype>& a, xvector<int> ia, 
        xmatrix<utype>& covar, utype& chisq, 
        void (*funcs)(utype, xvector<utype>&));

  void GCD(const xmatrix<int>& ma,const xmatrix<int>& mb,xmatrix<int>& mgcd); //CO20191201
  void GCD(const xmatrix<int>& ma,const xmatrix<int>& mb,xmatrix<int>& mgcd,xmatrix<int>& mx,xmatrix<int>& my); //CO20191201

  //AS20200512 BEGIN
  /// gives maximum among diagonal elements
  template<class utype> utype
    maxDiagonalElement(const xmatrix<utype>& M)
    {
      utype maxM = M[M.lrows][M.lcols];
      for (int i=M.lrows+1,j=M.lcols+1; i<=M.urows && j<=M.ucols; i++, j++){
        if (M[i][j] > maxM) maxM = M[i][j];
      }
      return maxM;
    }
  /// gives minimum among diagonal elements
  template<class utype> utype
    minDiagonalElement(const xmatrix<utype>& M)
    {
      utype minM = M[M.lrows][M.lcols];
      for (int i=M.lrows+1,j=M.lcols+1; i<=M.urows && j<=M.ucols; i++, j++){
        if (M[i][j] < minM) minM = M[i][j];
      }
      return minM;
    }
  //AS20200512 END

  //DX20210111 - added polarDecomposition() function
  template<class utype> void
    polarDecomposition(const xmatrix<utype>& transformation_matrix,
        xmatrix<utype>& rotation,
        xmatrix<utype>& deformation,
        bool check_orthogonal_rotation=true); __xprototype;

  template<class utype> void equilibrateMatrix(const xmatrix<utype>& A_orig, xmatrix<utype>& A, xmatrix<utype>& R, xmatrix<utype>& C, uint niter=100, utype tol=(utype)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_); __xprototype; //SD20220425  //CO20230313 - utype clang warnings

}

namespace aurostd {
  // ORTHOGONALITY
  template<class utype> utype
    orthogonality_defect(const xmatrix<utype>& basis) __xprototype;   // imported May 2018
  template<class utype> bool
    gaussian_reduce_two_vectors(xvector<utype>& B, xvector<utype>& C,utype eps) __xprototype;   // imported May 2018
  template<class utype> void
    reduce_A_in_ABC(xvector<utype>& A, xvector<utype>& B, xvector<utype>& C,utype eps) __xprototype;   // imported May 2018
  template<class utype> utype
    reduce_to_shortest_basis(const xmatrix<utype>& IN,xmatrix<utype>& OUT,utype eps,bool VERBOSE) __xprototype;   // imported May 2018
  template<class utype> xmatrix<utype>
    reduce_to_shortest_basis(const xmatrix<utype>& IN,utype eps,bool VERBOSE) __xprototype;   // imported May 2018
  template<class utype> xmatrix<utype>
    reduce_to_shortest_basis(const xmatrix<utype>& IN) __xprototype;   // imported May 2018
}

namespace aurostd {
  // EIGENPROBLEMS
  template<class utype>
    int jacobi(const xmatrix<utype> &ain,xvector<utype> &d,xmatrix<utype> &v) __xprototype;
  // Computes all eigenvalues and eigenvectors of a real symmetric xmatrix a[1..n][1..n].
  // On output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
  // v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
  // a. The function returns the number of Jacobi rotations that were required.
  //ME20190815
  template<class utype>
    xvector<utype> jacobiHermitian(xmatrix<xcomplex<utype> >&, char=0);
  template<class utype>
    xvector<utype> jacobiHermitian(xmatrix<xcomplex<utype> >&, xmatrix<xcomplex<utype> >&, char=0);
  template<class utype>
    void eigsrt(xvector<utype> &d,xmatrix<utype> &v) __xprototype;
  // Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi
  // or tqli,this routine sorts the eigenvalues into descending order, and rearranges
  // the columns of v correspondingly. The method is straight insertion and is N2 rather than NlogN;
  // but since you have just done an N3 procedure to get the eigenvalues, you can afford yourself
  // this little indulgence.
  template<class utype> 
    void QRDecomposition_HouseHolder(const xmatrix<utype>& mat,xmatrix<utype>& Q,xmatrix<utype>& R,utype tol=(utype)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_); //CO20191110 //CO20230313 - utype clang warnings
  template<class utype> 
    void QRDecomposition_HouseHolder_MW(const xmatrix<utype>& mat,xmatrix<utype>& Q,xmatrix<utype>& R,utype tol=(utype)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_); //CO20191110  //CO20230313 - utype clang warnings
  template<class utype> 
    void QRDecomposition_HouseHolder_TB(const xmatrix<utype>& mat,xmatrix<utype>& Q,xmatrix<utype>& R,utype tol=(utype)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_); //CO20191110  //CO20230313 - utype clang warnings
  // general Householder, mat is mxn, m>=n
  // See Numerical Linear Algebra, Trefethen and Bau, pg. 73
  template<class utype> 
    void getEHermite(utype a,utype b,xmatrix<utype>& ehermite); //CO+YL20191201
  template<class utype> 
    void getSmithNormalForm(const xmatrix<utype>& A,xmatrix<utype>& U,xmatrix<utype>& V,xmatrix<utype>& S,double tol=(double)_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_); //CO+YL20191201 //CO20230313 - utype clang warnings
  template<class utype>
    void tred2(const xmatrix<utype> &a,xvector<utype> &d,xvector<utype> &e) __xprototype;
  // Householder reduction of a real, symmetric matrix a[1..n][1..n].
  // On output, a is replaced by the orthogonal matrix Q eﬀecting the
  // transformation. d[1..n] returns the diagonal elments of
  // the tridiagonal matrix, and e[1..n] the oﬀ-diagonal elements, with e[1]=0.
  template<class utype>
    void tqli(xvector<utype> &d,xvector<utype> &e,xmatrix<utype> &z) __xprototype;
  // QL algorithm with implicit shifts, to determine the eigenvalues
  // and eigenvectors of a real, symmetric, tridiagonal matrix, or of a real,
  // symmetric matrix previously reduced by tred2
  // On input, d[1..n] contains the diagonal elements of the tridiagonal
  // matrix. On output, it returns the eigenvalues. The vectore[1..n]
  // inputs the subdiagonal elements of the tridiagonal matrix, with e[1] arbitrary.
  // On output e is destroyed. When finding only the eigenvalues, several lines
  // maybe omitted, as noted in the comments. If the eigenvectors of a tridiagonal
  // matrix are desired, the matrix z[1..n][1..n] is input as the identity
  // matrix. If the eigenvectors of a matrix that has been reduced by tred2
  // are required, then z is input as the matrix output by tred2.
  // In either case, the kth column of z returns the normalized eigenvector
  // corresponding to d[k].
  template<class utype>
    void balanc(xmatrix<utype> &a) __xprototype;
  // Given a matrix a[1..n][1..n], this routine replaces it by
  // a balanced matrix with i dentical eigenvalues. A symmetric matrix
  // is already balanced and is unaﬀected by this procedure. The
  // parameter RADIX should be the machine’s ﬂoating-point radix.
  template<class utype>
    void elmhes(xmatrix<utype> &a) __xprototype;
  // Reduction to Hessenberg form by the elimination method.
  // The real, nonsymmetric matrix a[1..n][1..n] is replaced by an upper
  // Hessenberg matrix with identical eigenvalues.
  // Recommended, but not required, is that this routine be preceded
  // by balanc. On output, the Hessenberg matrix is in elements a[i][j] with i<=j+1.
  // Elements with i>j+1 are to be thought of as zero,
  // but are returned with random values.
  template<class utype>
    void hqr(xmatrix<utype> &a,xvector<utype> &wr,xvector<utype> &wi) __xprototype;
  // Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n].
  // On input a can be exactly as output from elmhes; on output it is destroyed.
  // The real and imaginary parts of the eigenvalues are returned in
  //wr[1..n] and wi[1..n], respectively.
  template<class utype>
    void eigen(const xmatrix<utype> &ain,xvector<utype> &wr,xvector<utype> &wi) __xprototype;
  // Finds all eigenvalues of matrix a[1..n][1..n]. The real and imaginary parts
  // of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
  //ME20190718 BEGIN
  template<class utype>
    utype l1_norm(const xmatrix<utype>&);
  template<class utype>
    double l2_norm(const xmatrix<utype>&);
  template<class utype>
    double frobenius_norm(const xmatrix<utype>&);
  template<class utype>
    utype linf_norm(const xmatrix<utype>&);
  //ME20190718 END
}


// -------------------------------------------------------------- class cematrix
namespace aurostd {
  // namespace aurostd
  class cematrix {
    // add the least square and singular decomposition methods to xmatrix class
    private:
      int nrow, ncol; // number of rows and columes of matrix A
      xmatrix<double> M; // original matrix
      xvector<double> W; // A = U Z V^T
      xmatrix<double > U, V;
      xvector<double> a_vec; // a vector
      vector<double> a_nvec; // a vector (same as a_vec but in c++ vector form)
      double chisq; // value of the norm square
      xmatrix<double> Cov; // covariance matrix    
    public:
      cematrix();
      cematrix(const xmatrix<double> & A_in);
      ~cematrix();
      void LeastSquare(xvector<double> & y_vec, xvector<double> & y_sigma);
      void LeastSquare(xvector<double> & y_vec); //AS20200811
      //void SVDcmp_NR(); // Singular Value Decomposition
      bool SVDcmp_NR(); // Singular Value Decomposition
      void SVDsolve(xvector<double>& b_vec); // Solve vector a
      void SVDFit(xvector<double>& x_sigma, xvector<double>& y); // Least Square fitting using SVD
      void SVDvar(); // get covariance matrix
      xvector<double> GetFitVector() { return a_vec;}; // output a_vec
      xmatrix<double> InverseMatrix(); // inverse of a matrix, not necessary square
      xvector<double> EigenValues(); // calculate eigenvalues of a matrix, not necessary square
      double Pythag2(double a, double b); // (a^2+b^2)^(1/2) avoiding overflow and underflow
      double _sign(double a, double b) { return b > -1? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
      vector<double> AVec() { return a_nvec;};
      double ChiSQ() { return chisq;};
  };
}

// -------------------------------------------------------------- CMdet
namespace aurostd {
  template<class utype> utype
    CMdet(const xmatrix<utype>& B);
}

// -------------------------------------------------------------- Rotation matrices
namespace aurostd {
  template<class utype> xmatrix<utype>
    getRotationMatrix3D(const xvector<utype>& a,const xvector<utype>& b);
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
//CO20200404 - moving matrix() from pflow to aurostd because it is templated
//doesn't compile otherwise

namespace aurostd {
  // Dane Morgan style,adjusted by Stefano Curtarolo
  // Bringing definitions inside the template helps
  // constructing the right templates.
  template<class utype>
    class matrix {
      public:
        //NECESSARY PUBLIC CLASS METHODS - START
        //constructors - START
        matrix(void);
        matrix(const int m);
        matrix(const int m,const int n);
        matrix(const int m,const int n,const utype& inutype);
        matrix(const int m,const std::vector<utype>& inutypevec);
        matrix(const matrix& b);
        //constructors - STOP
        ~matrix(void);                         // destructor
        void clear();
        //NECESSARY PUBLIC CLASS METHODS - STOP
        // accessors
        void print(void);
        uint size(void) const {return (uint) mat.size();}
        matrix<utype> transpose(void) const;
        //  matrix<utype>::iterator begin();
        //  matrix<utype>::iterator end();
        // operator
        std::vector<utype>& operator[] (const int i) {assert(i>=0 && i<=(int) mat.size()); return mat[i];}
        const std::vector<utype>& operator[] (const int i) const {assert(i>=0 && i<=(int) mat.size()); return mat[i];};
        const matrix<utype>& operator=(const matrix<utype>& b);
        // mutators
        void push_back(const std::vector<utype>& inutypevec) {mat.push_back(inutypevec);}
        void pop_back(void) {mat.pop_back();}
        void vecvec2mat(const std::vector<std::vector<utype> >& inVV);
        void vec2mat(const std::vector<utype>& inV);
        //[CO20200404 - OBSOLETE]void clear(void) {mat.clear();}
        void insert(const int& id,const std::vector<utype>& inV) {mat.insert(mat.begin()+id,inV);}
        void erase(const int id);
        void erase_col(const int id);
      private:
        std::vector<std::vector<utype> > mat;
        void free();
        void copy(const matrix& b);
    };
  template <class utype> matrix<utype>  xmatrix2matrix(const xmatrix<utype>& );
  template <class utype> xmatrix<utype> matrix2xmatrix(const matrix<utype>& );
}


// ----------------------------------------------------------------------------

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

