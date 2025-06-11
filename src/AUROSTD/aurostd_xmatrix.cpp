// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2018
// fixed for g++ 4.5 on Mar11 2014
// compiles in g++6 and 7 2018

#ifndef _AUROSTD_XMATRIX_CPP_
#define _AUROSTD_XMATRIX_CPP_

//#define _AUROSTD_XMATRIX_DEFAULT_SIZE_ 3
//#define _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_ double(1.0e-6)
//#define _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_ double(1.0e-6)

//#ifndef _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_
//#define _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_ 1.0e-6 //DX20171025
//#endif

#ifndef XXEND
#define XXEND 1
#endif

#ifndef _xmatrix_epsilon
#define _xmatrix_epsilon 1.0e-15
#endif

#define _exponential_convergence 1.0e-18

#ifndef _AUROSTD_XCOMPLEX_H_
#include "aurostd_xcomplex.h"
#endif
#ifndef _AUROSTD_XSCALAR_H_
#include "aurostd_xscalar.h"
#endif
#ifndef _AUROSTD_XVECTOR_H_
#include "aurostd_xvector.h"
#endif
#ifndef _AUROSTD_XMATRIX_H_
#include "aurostd_xmatrix.h"
#endif
#ifndef _AUROSTD_XTENSOR6_H_
#include "aurostd_xtensor.h"
#endif

#ifndef __XOPTIMIZE
#define _XMATRIX_CHECK_BOUNDARIES_
#endif

// ----------------------------------------------------------------------------
// --------------------------------------------------------------- constructors
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // constructor
    xmatrix<utype>::xmatrix(int nrh,int nch,int nrl,int ncl) : msize(0) {
      lrows=std::min(nrl,nrh); //if(!nrh||!nch) lrows=0; this messes up convasp
      urows=std::max(nrl,nrh); //if(!nrh||!nch) urows=0; this messes up convasp
      lcols=std::min(ncl,nch); //if(!nrh||!nch) lcols=0; this messes up convasp
      ucols=std::max(ncl,nch); //if(!nrh||!nch) ucols=0; this messes up convasp
      refresh();  //CO20191112
#ifdef _XMATH_DEBUG_CONSTRUCTORS
      cerr << "M -> constructor:"
        << "  lrows=" << lrows << ", urows=" << urows
        << ", lcols=" << lcols << ", ucols=" << ucols
        << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
      if(msize>0) {
        corpus=new utype *[rows+XXEND]();  //HE20220613 initialize corpus memory
        if(!corpus){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::xmatrix<utype>::xmatrix():","allocation failure 1 (int,int,int,int)",_ALLOC_ERROR_);}
        corpus+= -lrows+ XXEND;
        corpus[lrows]= new utype[rows*cols+XXEND]();  //HE20220613 initialize corpus memory
        if(!corpus[lrows]){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::xmatrix<utype>::xmatrix():","allocation failure 2 (int,int,int,int)",_ALLOC_ERROR_);}
        corpus[lrows]+= -lcols+XXEND;
        for(int i=lrows+1;i<=urows;i++){corpus[i]=corpus[i-1]+cols;}  //this propagates previous line to all lrows
        clear();
      }
#ifdef _XMATH_DEBUG_CONSTRUCTORS
      cerr << "issquare=" << issquare << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex
        << ", sizeof=" << size << ", msize=" << msize << endl;
#endif
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                       // copy constructor
    xmatrix<utype>::xmatrix(const xmatrix<utype>& b) {
      init();
      copy(b);
    }  //CO20191112
  template<class utype>                                       // copy constructor
    xmatrix<utype>::xmatrix(const xvector<utype>& b) {
      init();
      copy(b);
    }  //CO20191112
  template<class utype>  // initializer_list constructor //HE20220616
    xmatrix<utype>::xmatrix(const std::initializer_list<std::initializer_list<utype>> ll) {
      // usage: xmatrix<double> new_matrix({{1,0, 2.0, 3.0, 4.0},
      //                                    {5,0, 6.0, 7.0, 8.0}});
      init();
      copy(ll);
    }

  template<class utype>
    void xmatrix<utype>::init(){  //HE20220613 initialize all members of xmatrix
      rows=0; lrows=0; urows= 0;
      cols=0; lcols=0; ucols=0;
      issquare=false; isfloat=false; iscomplex=false;
      size=0; msize=0;
    }

}

namespace aurostd {  // namespace aurostd
  template<class utype>                                       // copy constructor
    xmatrix<utype>::xmatrix(int vrows,int vcols, utype* a) : msize(0) {  // a starts from 0..
      lrows=1;urows=vrows;
      lcols=1;ucols=vcols;
      refresh();  //CO20191112
#ifdef _XMATH_DEBUG_CONSTRUCTORS
      cerr << "M -> copy constructor:"
        << "  lrows=" << lrows << ", urows=" << urows
        << ", lcols=" << lcols << ", ucols=" << ucols
        << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
      if(msize>0) {
        corpus=new utype *[rows+XXEND]();  //HE20220613 initialize corpus memory
        if(!corpus){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::xmatrix<utype>::xmatrix():","allocation failure 1 (int,int,utype*)",_ALLOC_ERROR_);}
        corpus+= -lrows+ XXEND;
        corpus[lrows]= new utype[rows*cols+XXEND]();  //HE20220613 initialize corpus memory
        if(!corpus[lrows]){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::xmatrix<utype>::xmatrix():","allocation failure 2 (int,int,utype*)",_ALLOC_ERROR_);}
        corpus[lrows]+= -lcols+XXEND;
        int i=0,j=0;
        for(i=lrows+1;i<=urows;i++){corpus[i]=corpus[i-1]+cols;}  //this propagates previous line to all lrows
        for(i=lrows;i<=urows;i++){
          for(j=lcols;j<=ucols;j++){
            corpus[i][j]=(utype) a[(i-1)*ucols+(j-1)]; // a.corpus[i][j];  // LIKE FORTRAN
          }
        }
        //      delete [] a;
      }
#ifdef _XMATH_DEBUG_CONSTRUCTORS
      cerr << "issquare=" << issquare << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex
        << ", sizeof=" << size << ", msize=" << msize << endl;
#endif
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------- destructor

namespace aurostd {  // namespace aurostd
  template<class utype>                                     // default destructor
    xmatrix<utype>::~xmatrix() {
      // cerr << "problem destructor xmatrix [1]" << endl;
      // free a xmatrix allocated with xmatrix()
#ifdef _XMATH_DEBUG_DESTRUCTORS
      cerr << "M -> default destructor:"
        << "  lrows=" << lrows << ", urows=" << urows
        << ", lcols=" << lcols << ", ucols=" << ucols
        << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
      free(); //CO20190808
    }
}
// ----------------------------------------------------------------------------
// -------------------------------------------------------- assigment operators

namespace aurostd { // namespace aurostd
  template<class utype>
    void xmatrix<utype>::free() { //CO20190808
      if(msize>0) {
        delete [] (corpus[lrows]+lcols-XXEND);
        delete [] (corpus+lrows-XXEND);
      }
      lrows=urows=lcols=ucols=0;refresh();
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
    void xmatrix<utype>::copy(const xmatrix<utype>& b) { //CO20190808
      if(lrows!=b.lrows||urows!=b.urows||lcols!=b.lrows||ucols!=b.ucols||msize!=b.msize) {    // if dims(this)!=dims(a) => build a new xmatrix !!!  //CO20190808 - VERY IMPORTANT that we not only check lrows/urows/lcols/ucols, but msize, as xmatrix could just have been initialized (msize==0)
        free();
        lrows=b.lrows;urows=b.urows;rows=b.rows;
        lcols=b.lcols;ucols=b.ucols;cols=b.cols;
        //[simply copy instead]refresh();
        issquare=bool(rows == cols);
        isfloat=_isfloat((utype) 0);
        iscomplex=_iscomplex((utype) 0);
        size=(char) sizeof(utype);
        msize=(long int) size*rows*cols;
#ifdef _XMATH_DEBUG_OPERATORS
        cerr << "M -> operator =::"
          << "  lrows=" << lrows << ", urows=" << urows
          << ", lcols=" << lcols << ", ucols=" << ucols
          << ", rows="  << rows  << ", cols="  << cols << endl;
#endif
        if(msize>0) {
          corpus=new utype *[rows+XXEND]();  //HE20220613 initialize corpus memory
          if(!corpus){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::xmatrix<utype>::copy():","allocation failure 1 in COPY",_ALLOC_ERROR_);}
          corpus+= -lrows+ XXEND;
          corpus[lrows]= new utype[rows*cols+XXEND]();  //HE20220613 initialize corpus memory
          if(!corpus[lrows]){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::xmatrix<utype>::copy():","allocation failure 2 in COPY",_ALLOC_ERROR_);}
          corpus[lrows]+= -lcols+XXEND;
          for(int i=lrows+1;i<=urows;i++){corpus[i]=corpus[i-1]+cols;}  //this propagates previous line to all lrows
        }
#ifdef _XMATH_DEBUG_CONSTRUCTORS
        cerr << "issquare=" << issquare << ", isfloat=" << isfloat << ", iscomplex=" << iscomplex
          << ", sizeof=" << size << ", msize=" << msize << endl;
#endif
      }
      if(corpus!=b.corpus){
        int i=0,j=0;
        for(i=0;i<rows;i++){
          for(j=0;j<cols;j++){
            corpus[i+lrows][j+lcols] = b.corpus[i+b.lrows][j+b.lcols];
          }
        }
      } //CO20190808 - we definitely have corpus now
    }
  template<class utype>
    void xmatrix<utype>::copy(const xvector<utype>& b) { //CO20190808
      xmatrix<utype> a(b.urows,1,b.lrows,1);
      for(int i=b.lrows;i<=b.urows;i++){a[i][1]=b[i];}
      copy(a);
    }
  template<class utype>
    void xmatrix<utype>::copy(std::initializer_list<std::initializer_list<utype>> ll) { //HE20220616
      int ll_rows = ll.size();
      int ll_cols = ll.begin()->size();
      xmatrix<utype> a(ll_rows,ll_cols,1,1);
      size_t new_row = 0;
      size_t new_col = 0;
      for (il2i l=ll.begin(); l<ll.end(); l++){
        new_row += 1;
        if (ll_cols != (int) l->size()) {
          stringstream message;
          message << "failure in copy - column size size mismatch ";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
        }
        for (ili entry = l->begin(); entry < l->end(); entry++) {
          new_col += 1;
          a[new_row][new_col] = *entry;
        }
        new_col = 0;
      }
      copy(a);
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                             // operator =
    xmatrix<utype>& xmatrix<utype>::operator=(const xmatrix<utype>& b) {  //CO20191112
      if(this!=&b) {copy(b);}
      return *this;
    }
  template<class utype>                                             // operator =
    xmatrix<utype>& xmatrix<utype>::operator=(const std::initializer_list<std::initializer_list<utype>> ll) {  //CO20191112
      // usage: xmatrix<double> new_matrix;
      // new_matrix = {{1,0, 2.0, 3.0, 4.0},
      //               {5,0, 6.0, 7.0, 8.0}};
      copy(ll);
      return *this;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
    void xmatrix<utype>::refresh(void) { //CO20190808
      rows=urows-lrows+1;      //if(!nrh||!nch) rows=0; this messes up convasp
      cols=ucols-lcols+1;      //if(!nrh||!nch) cols=0; this messes up convasp
      if(rows==0||cols==0) {
        cerr << "XMATRIX constructor: creating EMPTY xmatrix<utype>" << endl;
        lrows=0;urows=0;rows=0;
        lcols=0;ucols=0;cols=0;
      };
      issquare=bool(rows == cols);
      isfloat=_isfloat((utype) 0);
      iscomplex=_iscomplex((utype) 0);
      size=(char) (sizeof(utype));
      msize=(long int) size*rows*cols;
    }

  /// @brief shift the lower bounds of an xmatrix to new values
  /// \param new_lrows new lower bound for rows
  /// \param new_lcols new lower bound for columns
  ///
  /// @authors
  /// @mod{HE,20220915,created function}
  /// @note modifies the pointers but not the underlying data field
  template<class utype>
  void xmatrix<utype>::shift(int new_lrows, int new_lcols) {
    if (new_lrows == lrows && new_lcols == lcols) return;
    int row_shift = lrows - new_lrows;
    int col_shift = lcols - new_lcols;
    corpus += row_shift;
    lrows = new_lrows;
    lcols = new_lcols;
    urows = new_lrows + rows - 1;
    ucols = new_lcols + cols - 1;
    corpus[lrows]+= col_shift;
    for(int i=lrows+1;i<=urows;i++){corpus[i]=corpus[i-1]+cols;}
  }

  /// @brief shift the lower row bound of an xmatrix to a new value
  /// \param new_lrows new lower bound for rows
  /// @note modifies the pointers but not the underlying data field
  template<class utype>
  void xmatrix<utype>::shiftrow(int new_lrows) {
    shift(new_lrows, lcols);
  }

  /// @brief shift the lower column bound of an xmatrix to a new value
  /// \param new_lcols new lower bound for columns
  /// @note modifies the pointers but not the underlying data field
  template<class utype>
  void xmatrix<utype>::shiftcol(int new_lcols) {
    shift(lrows, new_lcols);
  }

  /// @brief reshape a xmatrix to given dimensions
  /// @param new_rows row count
  /// @param new_cols column count
  ///
  /// @authors
  /// @mod{HE,20220915,created function}
  /// @note When new_rows or new_cols is set to 0, it is automatically calculated.
  /// This is also true when new_cols is not given.
  /// @note lower bound for rows and columns is not changed
  /// @note modifies the pointers but not the underlying data field
  template<class utype>
  void xmatrix<utype>::reshape(uint new_rows, uint new_cols) {
    if (new_rows+new_cols == 0 ) {
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "_xmatrix<utype>::reshape() failed! rows or cols need to be defined", _INDEX_BOUNDS_);
    }
    if (new_rows == 0) new_rows = (rows*cols)/new_cols;
    if (new_cols == 0) new_cols = (rows*cols)/new_rows;
    if (new_rows < 1 || new_cols < 1) {
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "New dimensions cannot be less than one", _VALUE_ILLEGAL_);
    }
    int new_urows = lrows + new_rows - 1;
    int new_ucols = lcols + new_cols - 1;
    reshape(new_urows, new_ucols, lrows, lcols);
  }

  /// @brief reshape a xmatrix to given dimensions
  /// @param new_urows new upper bound for rows
  /// @param new_ucols new upper bound for columns
  /// @param new_lrows new lower bound for rows
  /// @param new_lcols new lower bound for columns
  ///
  /// @authors
  /// @mod{HE,20220915,created function}
  /// @note modifies the pointers but not the underlying data field
  template<class utype>
  void xmatrix<utype>::reshape(int new_urows, int new_ucols, int new_lrows, int new_lcols){
    int new_rows = new_urows - new_lrows + 1;
    int new_cols = new_ucols - new_lcols + 1;
    if (new_rows < 1 || new_cols < 1) {
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "New dimensions cannot be less than one", _VALUE_ILLEGAL_);
    }
    if (rows*cols != new_rows*new_cols) {
      stringstream message;
      message << "_xmatrix<utype>::reshape() failed! (original size: " << rows*cols << ") != (new size: " << new_rows*new_cols << ")";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    if ((lrows == new_lrows) && (lcols == new_lcols) && (urows == new_urows) && (ucols == new_ucols) ) return;
    // save the position of the main data
    utype* main_data_ptr =  corpus[lrows];
    // free the old row pointer storage
    delete [] (corpus+lrows-XXEND);

    // change the xmatrix properties
    urows = new_urows; ucols = new_ucols;
    lrows = new_lrows; lcols = new_lcols;
    refresh();

    // create new array to store row pointers
    corpus = new utype *[rows+XXEND]();

    // shift to fit lrows
    corpus += -lrows+ XXEND;
    // point the start back to the data
    corpus[lrows]= main_data_ptr;
    // fill the row pointer with start to the rows in the data array
    for(int i=lrows+1;i<=urows;i++){corpus[i]=corpus[i-1]+cols;}

  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------------ index operators
// ---------------------------------------------------------------- operator []

namespace aurostd {  // namespace aurostd
  template<class utype>
    // removed inline
    utype* xmatrix<utype>::operator[] (int ir) const {
#ifdef _XMATRIX_CHECK_BOUNDARIES_
      if(ir>urows)  {
        stringstream message;
        message << "_xmatrix<utype>_rows_high ir=" << ir << ", lrows=" << lrows << ", hrows=" << urows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      if(ir<lrows) {
        stringstream message;
        message << "_xmatrix<utype>_rows_low ir=" << ir << ", lrows=" << lrows << ", hrows=" << urows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
#endif
      return corpus[ir];
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator (i,j)
    // removed inline
    utype& xmatrix<utype>::operator()(int i,int j) const {
      //#ifndef XMATRIX_PERIODIC_BOUNDARY_CONDITIONS
#ifdef _XMATRIX_CHECK_BOUNDARIES_
      if(i>urows) {
        stringstream message;
        message << "M -> i=" << i << " > urows=" << urows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      if(i<lrows) {
        stringstream message;
        message << "M -> i=" << i << " < lrows=" << lrows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      if(j>ucols) {
        stringstream message;
        message << "M -> j=" << j << " > ucols=" << ucols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      if(j<lcols) {
        stringstream message;
        message << "M -> j=" << j << " < lcols=" << lcols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
#endif // _XMATRIX_CHECK_BOUNDARIES_
      return corpus[i][j];
    }
  ////#else
  //# ifdef XMATH_WARNING
  //# warning "XMATRIX_PERIODIC_BOUNDARY_CONDITIONS"
  //# endif
  //int ii=i,jj=j;
  //// if(ii>urows) ii=lrows+mod(i-lrows,urows-lrows+1);
  //// if(ii<lrows) ii=urows-mod(urows-i,urows-lrows+1);
  //// if(jj>ucols) jj=lcols+mod(j-lcols,ucols-lcols+1);
  //// if(jj<lcols) jj=ucols-mod(ucols-j,ucols-lcols+1);
  //if(ii>urows) ii-=rows;
  //if(ii<lrows) ii+=rows;
  //if(jj>ucols) jj-=cols;
  //if(jj<lcols) jj+=cols;
  //return corpus[ii][jj];
  //#endif
  //}
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                         // operator (i)
    xvector<utype> xmatrix<utype>::operator()(int i) const {
      xvector<utype> out(lcols,ucols);
      for(int j=lcols;j<=ucols;j++)
        out(j)=corpus[i][j];
      return out;
    }
}

//ME20180904 returns a matrix column as an xvector
namespace aurostd {  // namespace aurostd
  template<class utype>
    xvector<utype> xmatrix<utype>::getcol(int i) const {
      xvector<utype> out(lrows, urows);
      for (int j = lrows; j <= urows; j++) {
        out(j) = corpus[j][i];
      }
      return out;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
    xvector<utype> xmatrix<utype>::getdiag(int k,int _lrows) const { //CO20191210
      //first get length
      int _rows=0,i=0,j=0;
      for(i=lrows;i<=urows;i++){
        j=i+k;
        if(j>ucols){break;}
        _rows++;
      }
      xvector<utype> diag(_rows,_lrows);
      int index=_lrows;
      for(i=lrows;i<=urows;i++){
        j=i+k;
        if(j>ucols){break;}
        diag[index++]=corpus[i][j];
      }
      return diag;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype> void
    xmatrix<utype>::getxvecInPlace(xvector<utype>& xv_out, int lrow, int urow, int lcol, int ucol, int lrows_out) const {
      /// @brief Convert an xmatrix into an xvector given a set of indices. 
      ///
      /// @param xv_out the xvector that will be changed InPlace
      /// @param lrow lower row to include in xvector 
      /// @param urow upper row to include in xvector
      /// @param lcol lower column to include in xvector
      /// @param ucol upper column to include in xvector
      /// @param lrows_out the lower rows of the output xvector (if it does not match the xvector will be resized)
      ///
      /// @return void 
      ///
      /// @authors
      /// @mod{CO,2019110,created as getxmatInPlace (xvector overload) + xmatrix2xvector}
      /// @mod{AZ,20220711,refactored into getxvecInPlace}
      ///
      /// This is a function that slices an xmatrix into an xvector, It checks that the vector is 1-dimensional
      /// when slicing the xmatrix. Note that the indices are inclusive.
      /// @see
      /// @xlink{aurostd::getxmatInplace}
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      if(LDEBUG){
        cerr << "xmat=" << endl << (*this) << endl;
        cerr << "urows=" << urows << endl;
        cerr << "ucols=" << ucols << endl;
        cerr << "lrows=" << lrows << endl;
        cerr << "lcols=" << lcols << endl;
        cerr << "urow=" << urow << endl;
        cerr << "ucol=" << ucol << endl;
        cerr << "lrow=" << lrow << endl;
        cerr << "lcol=" << lcol << endl;
      }
      if(lrow<lrows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lrow<lrows",_INDEX_BOUNDS_);}
      if(urow>urows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"urow>urows",_INDEX_BOUNDS_);}
      if(lcol<lcols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lcol<lcols",_INDEX_BOUNDS_);}
      if(ucol>ucols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ucol>ucols",_INDEX_BOUNDS_);}
      if(lcol>ucol){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lcol>ucol",_INDEX_BOUNDS_);}
      if(lrow>urow){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lrow>urow",_INDEX_BOUNDS_);}

      if((ucol != lcol)&&(lrow != urow)){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"(ucol != lcol)&&(lrow != urow)",_INDEX_BOUNDS_);
      }

      int size_out = (ucol-lcol+1)*(urow-lrow+1);
      int urows_out = size_out+lrows_out-1;
      if(! (xv_out.rows==size_out && xv_out.lrows==lrows_out) ) { //check if necessary to create new object
        xvector<utype> xv(urows_out, lrows_out);
        xv_out=xv;
      }

      else if(ucol == lcol){
        for(int i = lrows_out; i <= urows_out; i++){
          xv_out(i) = corpus[lrow+i-1][lcol];
        }
        return;
      }
      for(int j = lrows_out; j <= urows_out; j++){
        xv_out(j) = corpus[lrow][lcol+j-1];
      }
    }
}
namespace aurostd {  // namespace aurostd
  /// @brief Convert an xmatrix into an xvector given a set of indices. 
  ///
  /// @param xv_out the xvector that will be changed InPlace
  /// @param lrow lower row to include in xvector 
  /// @param urow upper row to include in xvector
  /// @param lcol lower column to include in xvector
  /// @param ucol upper column to include in xvector
  ///
  /// @return xvector 
  ///
  /// This function allocates the xvector then calls getxvecInPlace(). 
  ///
  /// @authors
  /// @mod{CO,2019110,created as getxmatInPlace() + xmatrix2xvector()}
  /// @mod{AZ,20220711,refactored into getxvecInPlace()}
  ///
  /// @see
  /// @xlink{aurostd::getxvecInPlace()}
  /// @xlink{aurostd::getxmatInPlace()} 
  template<class utype> xvector<utype> 
    xmatrix<utype>::getxvec(int lrow, int urow, int lcol, int ucol, int lrows_out) const {
      xvector<utype> xv_out;
      (*this).getxvecInPlace(xv_out, lrow, urow, lcol, ucol, lrows_out);
      return xv_out;
    }
}
namespace aurostd {  // namespace aurostd
  /// @brief Convert xmatrix into xvector given a set of indices. 
  ///
  /// @param void 
  ///
  /// @return xvector 
  ///
  /// @authors
  /// @mod{AZ,20220711,created}
  ///
  /// This function performs the same task as getxvecInPlace(), but it 
  /// assumes that you have already passed it a 1-d slice of the xmatrix
  /// and then performs the type conversion.
  /// @see
  /// @xlink{aurostd::getxvecInPlace()}
  /// @xlink{aurostd::getxmatInPlace()} 
  template<class utype> xvector<utype>
    xmatrix<utype>::getxvec() const {
      return (*this).getxvec(lrows, urows, lcols, ucols);
    }
}
//CO20190808
namespace aurostd {  // namespace aurostd
  /// @brief Convert xmatrix into a submatrix given a set of indices. 
  ///
  /// @param mat_out the matrix that will be changed InPlace
  /// @param lrow lower row to include in submatrix 
  /// @param urow upper row to include in submatrix
  /// @param lcol lower column to include in submatrix
  /// @param ucol upper column to include in submatrix
  /// @param lrows_out the lower row of the matrix that will be written over mat_out (if it is not the same it will allocate a new matrix)
  /// @param lcols_out the lower column of the matrix that will be written over mat_out (if it is not the same it will allocate a new matrix)
  ///
  /// @return void 
  ///
  /// @authors
  /// @mod{CO,2019110,created}
  /// @mod{AZ,20220711,modified}
  ///
  /// This is a function that slices an xmatrix into into a submatrix, originally written
  /// by CO. When slicing the xmatrix. lrow is lower row. urow is upper row. lcol is 
  /// lower column and ucol is upper column. Note that the indices are inclusive.
  /// i.e any index given will be returned. Also lcol == ucol or lrow == urow in
  /// order to be a vector. This function is designed so the base "InPlace" function
  /// drives the rest of the functions. The InPlace designation means you must supply 
  /// a matrix that will then be changed InPlace to the matrix elements that are specified
  /// by the arguments.
  /// @note
  ///
  /// CO: this function returns a submatrix mat_out spanning urow:lrow,ucol:lcol of the original matrix
  /// lrows_out,lcols_out specifies lrows,lcols of mat_out
  /// this is different than submatrix(), which returns back a submatrix by cutting out irow,jcol
  template<class utype> void
    xmatrix<utype>::getxmatInPlace(xmatrix<utype>& mat_out,int lrow,int urow,int lcol,int ucol,int lrows_out,int lcols_out) const { //lrow, lcol references corpus, lrows_out references output  //CO20191110
      //AZ20220627 START
      if(lrow<lrows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lrow<lrows",_INDEX_BOUNDS_);}
      if(urow>urows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"urow>urows",_INDEX_BOUNDS_);}
      if(lcol<lcols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lcol<lcols",_INDEX_BOUNDS_);}
      if(ucol>ucols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ucol>ucols",_INDEX_BOUNDS_);}
      if(lcol>ucol){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lcol>ucol",_INDEX_BOUNDS_);}
      if(lrow>urow){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lrow>urow",_INDEX_BOUNDS_);}
      //AZ20220627 END
      int rows_out=(urow-lrow)+1;
      int cols_out=(ucol-lcol)+1;

      if(! ( mat_out.rows==rows_out && mat_out.cols==cols_out && mat_out.lrows==lrows_out && mat_out.lcols==lcols_out ) ) { //check if necessary to create new object
        xmatrix<utype> mat(rows_out,cols_out,lrows_out,lcols_out);
        mat_out=mat;
      }
      for(int i=lrow;i<=urow;i++){
        for(int j=lcol;j<=ucol;j++){mat_out[(i-lrow)+mat_out.lrows][(j-lcol)+mat_out.lcols]=corpus[i][j];}
      }
    }
  template<class utype> xmatrix<utype>
    xmatrix<utype>::getxmat(int lrow,int urow,int lcol,int ucol,int lrows_out,int lcols_out) const { //lrow, lcol references corpus, lrows_out references output  //CO20191110
      xmatrix<utype> xmat;
      (*this).getxmatInPlace(xmat,lrow,urow,lcol,ucol,lrows_out,lcols_out);
      return xmat;
    }
}

//CO20190808
namespace aurostd {  // namespace aurostd
  template<class utype>
    void xmatrix<utype>::setrow(const xvector<utype>& row,int irow) {  //CO20191110
      return setmat(row,irow,false);
      //[OVERLOAD WITH SETMAT()]if(row.lrows!=lcols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"row.lrows!=lcols",_INPUT_ILLEGAL_);}
      //[OVERLOAD WITH SETMAT()]if(row.urows!=ucols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"row.urows!=ucols",_INPUT_ILLEGAL_);}
      //[OVERLOAD WITH SETMAT()]if(irow<lrows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"irow<lrows",_INPUT_ILLEGAL_);}
      //[OVERLOAD WITH SETMAT()]if(irow>urows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"irow>urows",_INPUT_ILLEGAL_);}
      //[OVERLOAD WITH SETMAT()]for(int j=row.lrows;j<=row.urows;j++){corpus[irow][j]=row[j];}
    }
}

//CO20190808
namespace aurostd {  // namespace aurostd
  template<class utype>
    void xmatrix<utype>::setcol(const xvector<utype>& col,int icol) {  //CO20191110
      return setmat(col,icol,true);
      //[OVERLOAD WITH SETMAT()]if(col.lrows!=lrows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"col.lrows!=lrows",_INPUT_ILLEGAL_);}
      //[OVERLOAD WITH SETMAT()]if(col.urows!=urows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"col.urows!=urows",_INPUT_ILLEGAL_);}
      //[OVERLOAD WITH SETMAT()]if(icol<lcols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"icol<lcols",_INPUT_ILLEGAL_);}
      //[OVERLOAD WITH SETMAT()]if(icol>ucols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"icol>ucols",_INPUT_ILLEGAL_);}
      //[OVERLOAD WITH SETMAT()]for(int j=col.lrows;j<=col.urows;j++){corpus[j][icol]=col[j];}
    }
}

//CO20190808
namespace aurostd {  // namespace aurostd
  template<class utype>
    void xmatrix<utype>::setmat(const xmatrix<utype>& mat,int lrow,int lcol) { //these are the starting lrow, lcol, end is dictated by size of mat //CO20191110
#ifdef _XMATRIX_CHECK_BOUNDARIES_
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      int urow=lrow+mat.rows-1; //ending row
      int ucol=lcol+mat.cols-1; //ending col
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " urow=" << urow << endl;
        cerr << __AFLOW_FUNC__ << " ucol=" << ucol << endl;
      }
      if(lrow<lrows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lrow<lrows",_VALUE_ILLEGAL_);}
      if(urow>urows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"urow>urows",_VALUE_ILLEGAL_);}
      if(lcol<lcols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lcol<lcols",_VALUE_ILLEGAL_);}
      if(ucol>ucols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ucol>ucols",_VALUE_ILLEGAL_);}
#endif
      for(int i=mat.lrows;i<=mat.urows;i++){
        for(int j=mat.lcols;j<=mat.ucols;j++){corpus[lrow+i-mat.lrows][lcol+j-mat.lcols]=mat[i][j];}
      }
    }
  template<class utype>
    void xmatrix<utype>::setmat(const xvector<utype>& xv,int icol,bool col) { //replace icol (col==true) or row (col==false) //CO20191110
      int lrow=1,lcol=1;
      if(col==true){
        lrow=lrows; //starting row
        lcol=icol; //starting col
      }else{
        lrow=icol; //starting row
        lcol=lcols; //starting col
      }
#ifdef _XMATRIX_CHECK_BOUNDARIES_
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      int urow=1,ucol=1;
      if(col==true){
        urow=lrow+xv.rows-1; //ending row
        ucol=icol; //ending col
      }else{
        urow=icol; //ending row
        ucol=lcols+xv.rows-1; //ending col
      }
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " lrow=" << lrow << endl;
        cerr << __AFLOW_FUNC__ << " urow=" << urow << endl;
        cerr << __AFLOW_FUNC__ << " lcol=" << lcol << endl;
        cerr << __AFLOW_FUNC__ << " ucol=" << ucol << endl;
      }
      if(lrow<lrows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lrow<lrows",_VALUE_ILLEGAL_);}
      if(urow>urows){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"urow>urows",_VALUE_ILLEGAL_);}
      if(lcol<lcols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"lcol<lcols",_VALUE_ILLEGAL_);}
      if(ucol>ucols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ucol>ucols",_VALUE_ILLEGAL_);}
#endif
      if(col==true){
        for(int i=xv.lrows;i<=xv.urows;i++){corpus[lrow+i-xv.lrows][icol]=xv[i];}
      }else{
        for(int i=xv.lrows;i<=xv.urows;i++){corpus[icol][lcol+i-xv.lrows]=xv[i];}
      }
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------- index operators with boundary conditions

namespace aurostd {  // namespace aurostd
  template<class utype>                        // operator () boundary conditions
    // removed inline
    utype& xmatrix<utype>::operator()(int i,int j,bool bc) const {
      if(bc==BOUNDARY_CONDITIONS_PERIODIC) {
        int ii=i,jj=j;
        if(ii==urows+1) ii=lrows; // fast switching
        if(ii==lrows-1) ii=urows; // fast switching
        if(ii>urows) ii=lrows+mod(i-lrows,urows-lrows+1);
        if(ii<lrows) ii=urows-mod(urows-i,urows-lrows+1);
        if(jj==ucols+1) jj=lcols; // fast switching
        if(jj==lcols-1) jj=ucols; // fast switching
        if(jj>ucols) jj=lcols+mod(j-lcols,ucols-lcols+1);
        if(jj<lcols) jj=ucols-mod(ucols-j,ucols-lcols+1);
#ifdef _XMATRIX_CHECK_BOUNDARIES_
        if(ii>urows) {
          stringstream message;
          message << "V -> ii=" << ii << " > urows" << urows << " <<  BC=" << bc;
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        if(ii<lrows) {
          stringstream message;
          message << "V -> ii=" << ii << " < lrows" << lrows << " <<  BC=" << bc;
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        if(jj>ucols) {
          stringstream message;
          message << "V -> jj=" << jj << " > ucols" << ucols << " <<  BC=" << bc;
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        if(jj<lcols) {
          stringstream message;
          message << "V -> jj=" << jj << " < lcols" << lcols << " <<  BC=" << bc;
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
#endif
        return corpus[ii][jj];
      }
      else { // ensure that this function always hit a return //HE20220616
#ifdef _XMATRIX_CHECK_BOUNDARIES_
        if(i>urows) {
          stringstream message;
          message << "M -> i=" << i << " > urows=" << urows;
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        if(i<lrows) {
          stringstream message;
          message << "M -> i=" << i << " < lrows=" << lrows;
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        if(j>ucols) {
          stringstream message;
          message << "M -> j=" << j << " > ucols=" << ucols;
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        if(j<lcols) {
          stringstream message;
          message << "M -> j=" << j << " < lcols=" << lcols;
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
#endif
        return corpus[i][j];
      }
    }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------- math unary operators

// -------------------------------------------------- operator xmatrix += xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>&
    // removed inline
    xmatrix<utype>::operator +=(const xmatrix<utype>& r)
    {
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator +=: ");
      printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
      printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
      printf("                 ");
      printf("r.lrows=%i, r.urows=%i, ",r.lrows,r.urows);
      printf("r.lcols=%i, r.ucols=%i\n",r.lcols,r.ucols);
#endif
      if(this->rows!=r.rows||this->cols!=r.cols) {
        string message = "(this->rows!=r.rows||this->cols!=r.cols)";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
          corpus[i+lrows][j+lcols]+=r[i+r.lrows][j+r.lcols];
      return *this;
    }

  template<class utype> xmatrix<utype>&
    xmatrix<utype>::operator +=(const std::initializer_list<std::initializer_list<utype>> ll){ //HE20220616
      int ll_rows = ll.size();
      int ll_cols = ll.begin()->size();
      if(this->rows!=ll_rows||this->cols!=ll_cols) {
        string message = "shape miss-match";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      size_t new_row = lrows;
      size_t new_col = lcols;
      for (il2i l=ll.begin(); l<ll.end(); l++){
        if (ll_cols != (int) l->size()) {
          string message = "column size size mismatch ";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
        }
        for (ili entry = l->begin(); entry < l->end(); entry++) {
          corpus[new_row][new_col] += *entry;
          new_col += 1;
        }
        new_row += 1;
      }
      return *this;
    }
}

// -------------------------------------------------- operator xmatrix -= xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>&
    // removed inline
    xmatrix<utype>::operator -=(const xmatrix<utype>& r)
    {
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator -=: ");
      printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
      printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
      printf("                 ");
      printf("r.lrows=%i, r.urows=%i, ",r.lrows,r.urows);
      printf("r.lcols=%i, r.ucols=%i\n",r.lcols,r.ucols);
#endif
      if(this->rows!=r.rows||this->cols!=r.cols) {
        string message = "(this->rows!=r.rows||this->cols!=r.cols)";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
          corpus[i+lrows][j+lcols]-=r[i+r.lrows][j+r.lcols];
      return *this;
    }
  template<class utype> xmatrix<utype>&
    xmatrix<utype>::operator -=(const std::initializer_list<std::initializer_list<utype>> ll) {//HE20220616
      int ll_rows = ll.size();
      int ll_cols = ll.begin()->size();
      if(this->rows!=ll_rows||this->cols!=ll_cols) {
        string message = "shape miss-match";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      size_t new_row = lrows;
      size_t new_col = lcols;
      for (il2i l=ll.begin(); l<ll.end(); l++){
        if (ll_cols != (int) l->size()) {
          string message = "column size size mismatch ";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
        }
        for (ili entry = l->begin(); entry < l->end(); entry++) {
          corpus[new_row][new_col] -= *entry;
          new_col += 1;
        }
        new_row += 1;
      }
      return *this;
    }

}

// -------------------------------------------------- operator xmatrix *= xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>&
    // removed inline
    xmatrix<utype>::operator *=(const xmatrix<utype>& b)
    {
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator *=: ");
      printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
      printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
      printf("                 ");
      printf("b.lrows=%i, b.urows=%i, ",b.lrows,b.urows);
      printf("b.lcols=%i, b.ucols=%i\n",b.lcols,b.ucols);
#endif
      if(!this->issquare||!b.issquare||this->rows!=b.rows)
        throw aurostd::xerror(__AFLOW_FILE__,"xmatrix<utype>::operator *=():","failure in operator*=: defined only for square xmatrixes with equal dimensions",_INPUT_ILLEGAL_);  //CO20191112

      xmatrix<utype> a(this->urows,this->ucols,this->lrows,this->lcols);
      int i=0,j=0,k=0,ii=0,jj=0,kk=0;
      utype *bk,*ai,aik,*thisi;

      for(i=this->lrows;i<=this->urows;i++)
        for(j=this->lcols;j<=this->ucols;j++) {
          a.corpus[i][j]=this->corpus[i][j];
          this->corpus[i][j]=(utype) 0;
        }
      for(i=this->lrows,ii=a.lrows;i<=this->urows;i++,ii++) {
        thisi=this->corpus[i];
        ai=a[ii];
        for(k=a.lcols,kk=b.lcols;k<=a.ucols;k++,kk++) {
          bk=b[kk];
          aik=ai[k];
          for(j=this->lrows,jj=b.lcols;j<=this->urows;j++,jj++)
            thisi[j]+=aik*bk[jj];
        }
      }
      return *this;
    }
}

// -------------------------------------------------- operator xmatrix *= utype
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>&
    // removed inline
    xmatrix<utype>::operator *=(utype r)
    {  //CO20191110
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator *=: ");
      printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
      printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
      printf("                 ");
      printf("b.lrows=%i, b.urows=%i, ",b.lrows,b.urows);
      printf("b.lcols=%i, b.ucols=%i\n",b.lcols,b.ucols);
#endif
      for(int i=lrows;i<=urows;i++)
        for(int j=lcols;j<=ucols;j++) {
          corpus[i][j]*=r;
        }
      return *this;
    }
}

// -------------------------------------------------- operator xmatrix /= utype
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>&
    // removed inline
    xmatrix<utype>::operator /=(utype r){  //CO20191110
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator *=: ");
      printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
      printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
      printf("                 ");
      printf("b.lrows=%i, b.urows=%i, ",b.lrows,b.urows);
      printf("b.lcols=%i, b.ucols=%i\n",b.lcols,b.ucols);
#endif
      for(int i=lrows;i<=urows;i++)
        for(int j=lcols;j<=ucols;j++) {
          corpus[i][j]/=r;
        }
      return *this;
    }
  template<class utype> xmatrix<utype>&
    // removed inline
    xmatrix<utype>::operator /=(const xmatrix<utype>& a){  //CO20191201 - right matrix division
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator *=: ");
      printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
      printf("this->lcols=%i, this->ucols=%i\n",this->lcols,this->ucols);
      printf("                 ");
      printf("b.lrows=%i, b.urows=%i, ",b.lrows,b.urows);
      printf("b.lcols=%i, b.ucols=%i\n",b.lcols,b.ucols);
#endif
      *this=*this*inverse(a);
      return *this;
    }
}

// ----------------------------------------------------------- operator +xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>
    xmatrix<utype> operator+(const xmatrix<utype>& a) {
      return a;
    }
}

// ----------------------------------------------------------- operator -xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>
    xmatrix<utype> operator-(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for (int i=a.lrows;i<=a.urows;i++)
        for (int j=a.lcols;j<=a.ucols;j++)
          c[i][j]=-a[i][j];
      return c;
    }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------ math binary operators

// ----------------------------------------------------------------------------
// --------------------------------------------------- operator xmatrix + xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>
    xmatrix<utype> operator+(const xmatrix<utype>& a,const xmatrix<utype>& b) {

#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator +: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
      printf("M -> operator +: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",b.lrows,b.urows,b.lcols,b.ucols);
#endif
      if(a.rows!=b.rows||a.cols!=b.cols) {
        string message = "(a.rows!=b.rows||a.cols!=b.cols)";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      int lr=1, ur=1, lc=1, uc=1;
      if ((a.lrows == b.lrows) && (a.lcols==b.lcols)) {lr = a.lrows; ur = a.urows; lc = a.lcols; uc = a.ucols;}
      else { ur = a.rows; uc = a.cols;}
      xmatrix<utype> c(ur, uc, lr, lc);
      int i,j;
      utype *bi,*ci,*ai;
      for(i=0;i<a.rows;i++) {
        ai=a[i+a.lrows];bi=b[i+b.lrows];ci=c[i+c.lrows];
        for(j=0;j<a.cols;j++)
          ci[j+c.lcols]=ai[j+a.lcols]+bi[j+b.lcols];}
      return c;
    }
}

// ----------------------------------------------------------------------------
// --------------------------------------------------- operator xmatrix - xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>
    xmatrix<utype> operator-(const xmatrix<utype>& a,const xmatrix<utype>& b) {

#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator +: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
      printf("M -> operator +: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",b.lrows,b.urows,b.lcols,b.ucols);
#endif
      if(a.rows!=b.rows||a.cols!=b.cols) {
        string message = "(a.rows!=b.rows||a.cols!=b.cols)";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      xmatrix<utype> c(a.rows,a.cols);
      int i,j;
      utype *bi,*ci,*ai;
      for(i=0;i<a.rows;i++) {
        ai=a[i+a.lrows];bi=b[i+b.lrows];ci=c[i+c.lrows];
        for(j=0;j<a.cols;j++)
          ci[j+c.lcols]=ai[j+a.lcols]-bi[j+b.lcols];};
      return c;
    }
}

// ----------------------------------------------------------------------------
// --------------------------------------------------- operator xmatrix * xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>
    xmatrix<utype> operator*(const xmatrix<utype>& a,const xmatrix<utype>& b) {
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator *: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
      printf("M -> operator *: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",b.lrows,b.urows,b.lcols,b.ucols);
#endif
      if(a.cols!=b.rows) {
        //ME20190814 - eliminate exit
        string message = "a.cols != b.rows";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      xmatrix<utype> c(a.rows,b.cols);
      int i=0,j=0,k=0,ii=0,jj=0,kk=0;
      // register
      utype *bk,*ci,*ai,aik;
      for(i=c.lrows,ii=a.lrows;i<=c.urows;i++,ii++) {
        ci=c[i];
        ai=a[ii];
        //for(k=a.lcols,kk=b.lcols;k<=a.ucols;k++,kk++)
        for(k=a.lcols,kk=b.lrows;k<=a.ucols;k++,kk++)
        { //CO20200106 - patching for auto-indenting
          bk=b[kk];
          aik=ai[k];
          for(j=c.lcols,jj=b.lcols;j<=c.ucols;j++,jj++)
            ci[j]+=aik*bk[jj];
        }
      }
      //for(i=c.lrows,ii=a.lrows;i<=c.urows;i++,ii++)          // 48% slower than the
      //for(k=a.lcols,kk=b.lrows;k<=a.ucols;k++,kk++)        // previous optimized
      //for(j=c.lcols,jj=b.lcols;j<=c.ucols;j++,jj++)      // routine
      //c[i][j]+=a[ii][k]*b[kk][jj];	
      //for(i=c.lrows,ii=a.lrows;i<=c.urows;i++,ii++)          // 66% slower than the
      //for(k=a.lcols,kk=b.lrows;k<=a.ucols;k++,kk++)        // previous optimized
      //for(j=c.lcols,jj=b.lcols;j<=c.ucols;j++,jj++)      // routine
      //c(i,j)+=a(ii,k)*b(kk,jj);
      return  c;
    }

  //ME20190814 - multiplication of a real matrix with a complex matrix
  template<class utype>
    xmatrix<xcomplex<utype> > operator*(const xmatrix<utype>& a, const xmatrix<xcomplex<utype> >& b) {
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator *: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
      printf("M -> operator *: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",b.lrows,b.urows,b.lcols,b.ucols);
#endif
      if (a.cols!=b.rows) {
        string message = "a.cols != b.rows";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }

      xmatrix<xcomplex<utype> > c(a.rows, b.cols);
      int i=0, j=0, k=0, ii=0, jj=0, kk=0;
      utype *ai, aik = (utype)0;
      xcomplex<utype> *bk, *ci;
      for (i = c.lrows, ii = a.lrows; i <= c.urows; i++, ii++) {
        ci = c[i];
        ai = a[ii];
        for (k = a.lcols, kk = b.lrows; k <= a.ucols; k++, kk++) {
          bk = b[kk];
          aik = ai[k];
          for (j = c.lcols, jj = b.lcols; j <= c.ucols; j++, jj++) {
            ci[j].re += aik * bk[jj].re;
            ci[j].im += aik * bk[jj].im;
          }
        }
      }
      return c;
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                               // operator xmatrix * xvector
    xvector<utype> operator*(const xmatrix<utype>& a,const xvector<utype>& b) {
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator *: a.lrows=%i, a.urows=%i, a.lcols=%i, a.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
      printf("M -> operator *: b.lrows=%i, b.urows=%i \n",b.lrows,b.urows);
#endif
      if(a.cols!=b.rows) {
        stringstream message;
        message << "xmatrix * xvector: Matrix and vector have different dimensions.";
        message << " a.cols = " << a.cols << ", b.rows = " << b.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      xvector<utype> c(a.lrows,a.urows);
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          c(i)+=a(i,j)*b(j-a.lcols+b.lrows); //HE20220912 xmatrix and xvector can start at different lcols/lrows
      return  c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                               // operator xvector * xmatrix
    xvector<utype> operator*(const xvector<utype>& a,const xmatrix<utype>& b) {
#ifdef _XMATH_DEBUG_OPERATORS
      printf("M -> operator *: a.lrows=%i, a.urows=%i \n",a.lrows,a.urows);
      printf("M -> operator *: b.lrows=%i, b.urows=%i, b.lcols=%i, b.ucols=%i\n",a.lrows,a.urows,a.lcols,a.ucols);
#endif
      if(a.rows!=b.rows) {
        stringstream message;
        message << "xvector * xmatrix: Vector and matrix have different dimensions.";
        message << " a.rows = " << a.rows << ", b.rows = " << b.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      xvector<utype> c(b.lcols,b.ucols);
      for(int i=b.lcols;i<=b.ucols;i++)
        for(int j=a.lrows;j<=a.urows;j++)
          //      c[i]+=a[j]*b[j-a.lrows+b.lrows][i];
          c(i)+=a(j)*b(j-a.lrows+b.lrows,i);
      return  c;
    }
}

//ME20200330 - Multiplication of a real matrix with a complex vector
namespace aurostd {
  template<class utype>
    xvector<xcomplex<utype> > operator*(const xmatrix<utype>& a, const xvector<xcomplex<utype> >& b) {
      if (a.cols != b.rows) {
        stringstream message;
        message << "xmatrix * xvector: Matrix and vector have different dimensions.";
        message << " a.cols = " << a.cols << ", b.rows = " << b.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      xvector<xcomplex<utype> > c(a.lrows, a.urows);
      for (int i = a.lrows; i <= a.urows; i++) {
        for (int j = a.lcols; j <= a.ucols; j++) {
          c[i].re += a[i][j] * b[j - a.lcols+b.lrows].re; //HE20220912 xmatrix and xvector can start at different lcols/lrows
          c[i].im += a[i][j] * b[j - a.lcols+b.lrows].im;
        }
      }
      return c;
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
//SD20220705
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                 // operator xmatrix + scalar
    operator+(const utype s,const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=(utype) a[i][j]+(utype) s;
      return c;
    }
}

// ----------------------------------------------------------------------------
//SD20220705
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                //  operator scalar + xmatrix
    operator+(const xmatrix<utype>& a,const utype s) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=(utype) s+(utype) a[i][j];
      return c;
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                 // operator xmatrix * scalar
    operator*(const utype s,const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=(utype) a[i][j]*(utype) s;
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                //  operator scalar * xmatrix
    operator*(const xmatrix<utype>& a,const utype s) {
      return s*a;
    }
}

//ME20200329 - real * complex matrix
namespace aurostd {
  template<class utype> xmatrix<xcomplex<utype> >
    operator*(utype s, const xmatrix<xcomplex<utype> >& a) {
      xmatrix<xcomplex<utype> > c(a.urows, a.ucols, a.lrows, a.lcols);
      for (int i = c.lrows; i <= c.urows; i++) {
        for (int j = c.lcols; j <= c.ucols; j++) {
          c[i][j].re = a[i][j].re * s;
          c[i][j].im = a[i][j].im * s;
        }
      }
      return c;
    }

  template<class utype> xmatrix<xcomplex<utype> >
    operator*(const xmatrix<xcomplex<utype> >& a, utype s) {
      return s*a;
    }

  template<class utype> xmatrix<xcomplex<utype> >
    operator/(const xmatrix<xcomplex<utype> >& a, utype s) {
      return ((utype) (1/s)) * a;
    }
}


// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                 // operator xmatrix / scalar
    operator/(const xmatrix<utype>& a,const utype s) {
      return (utype) ((utype)1/s)*a;                     //DX20170115 - add utype to 1/s to account for xcomplex
    }
  template<class utype> xmatrix<utype>                 // operator xmatrix / xmatrix
    operator/(const xmatrix<utype>& a,const xmatrix<utype>& b) {  //CO20191201
      return a*inverse(b);
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// CONDITIONALS

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    __identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_,const char& _mode_) {
      //ME20200725 - Changed exit to return false.
      if(a.rows!=b.rows) return false;
      if(a.cols!=b.cols) return false;
      bool output=TRUE;
      if(a.isfloat || a.iscomplex) {
        if(_mode_==1) { // relative tolerance
          for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
            for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++) {
              output=output*(((abs(a[i][j]-b[ii][jj]))/(abs(a[i][j])+abs(b[ii][jj])+_tol_))<=_tol_);
              if(output==FALSE) return (bool) output;
            }
        }
        if(_mode_==0) { // absolute tolerance (faster)  DEFAULT
          for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
            for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++) {
              output=output*(abs(a[i][j]-b[ii][jj])<=_tol_);
              if(output==FALSE) return (bool) output;
            }
        }
      } else {
        for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
          for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++) {
            output=output*(a[i][j]==b[ii][jj]);
            if(output==FALSE) return (bool) output;
          }
      }
      return (bool) output;
    }

  template<class utype> bool                             // is xmatrix == xmatrix ?
    identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_,const char& _mode_) {
      //ME20190814 BEGIN
      if ((a.rows != b.rows) || (a.cols != b.cols)) return false;
      //if(a.isfloat || a.iscomplex) { ME20190814 - this doesn't work for xcomplex because abs() and _tol_ have different types //[CO20200106 - close bracket for indenting]}
      //ME20190814 END
      if(a.isfloat) {
        if(_mode_==1) { // relative tolerance
          for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
            for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++)
              if((abs(a[i][j]-b[ii][jj])/(abs(a[i][j])/2.0+abs(b[ii][jj])/2.0+_tol_))>=_tol_) return FALSE;
        }
        if(_mode_==0) { // absolute tolerance (faster)  DEFAULT
          for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
            for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++)
              if(abs(a[i][j]-b[ii][jj])>=_tol_) return FALSE;
        }
      } else {
        for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
          for(int j=a.lcols,jj=b.lcols;j<=a.ucols;j++,jj++)
            if((a[i][j]!=b[ii][jj])) return FALSE;
      }
      return TRUE; // if FALSE has never found....
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_,(char) 0);  // relative
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    rel_identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_,(char) 1);  // relative
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    abs_identical(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_,(char) 0);  // absolute
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    identical(const xmatrix<utype>& a,const xmatrix<utype>& b) {
      return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    operator==(const xmatrix<utype>& a,const xmatrix<utype>& b) {
      return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
    isdifferent(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
      return (bool) !identical(a,b,_tol_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
    isdifferent(const xmatrix<utype>& a,const xmatrix<utype>& b) {
      return (bool) !identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    isequal(const xmatrix<utype>& a,const xmatrix<utype>& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    isequal(const xmatrix<utype>& a,const xmatrix<utype>& b) {
      return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
    operator!=(const xmatrix<utype>& a,const xmatrix<utype>& b) {
      return (bool) !identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  //ME20190814 - xcomplex version
  // namespace aurostd
  template<class utype> bool
    identical(const xmatrix<xcomplex<utype> >& a, const xmatrix<xcomplex<utype> >& b, const utype& _tol_, const char& _mode_) {
      if ((a.rows != b.rows) || (a.cols != b.cols)) return false;
      if (_mode_ == 1) {  // relative tolerance
        for (int i = a.lrows, ii = b.lrows; i <= a.urows; i++, ii++) {
          for (int j = a.lcols, jj = b.lcols; i <= a.ucols; j++, jj++) {
            if ((abs(a[i][j].re - b[ii][jj].re)/abs(a[i][j].re/2.0 + abs(b[ii][jj].re)/2.0 + _tol_)) >= _tol_) return false;
            if ((abs(a[i][j].im - b[ii][jj].im)/abs(a[i][j].im/2.0 + abs(b[ii][jj].im)/2.0 + _tol_)) >= _tol_) return false;
          }
        }
      } else if (_mode_ == 0) {  // absolute tolerance (faster) DEFAULT
        for (int i = a.lrows, ii = b.lrows; i <= a.urows; i++, ii++) {
          for (int j = a.lcols, jj = b.lcols; j <= a.ucols; j++, jj++) {
            if (isdifferent(a, b, _tol_)) return false;
          }
        }
      } else {  // unknown mode
        string message = "Unknown mode " + utype2string<char>(_mode_) + ".";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      return true;
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    identical(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_,(char) 0);  // relative
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    rel_identical(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_,(char) 1);  // relative
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    abs_identical(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_,(char) 0);  // absolute
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    identical(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b) {
      return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    operator==(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b) {
      return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
    isdifferent(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b,const utype& _tol_) {
      return (bool) !identical(a,b,_tol_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
    isdifferent(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b) {
      return (bool) !identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    isequal(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix == xmatrix ?
    isequal(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b) {
      return (bool) identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool                             // is xmatrix != xmatrix ?
    operator!=(const xmatrix<xcomplex<utype> >& a,const xmatrix<xcomplex<utype> >& b) {
      return (bool) !identical(a,b,(utype) _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_,(char) 0);
    }

  // namespace aurostd
  template<class utype> bool
    isinteger(const xmatrix<utype>& a,const utype& tol) {
      if(a.isfloat || a.iscomplex) {
        for(int i=a.lrows;i<=a.urows;i++)
          for(int j=a.lcols;j<=a.ucols;j++)
            if(isinteger(a[i][j],tol)==FALSE) return FALSE;
      }
      return TRUE;
    }

  // namespace aurostd
  //CO START
  template<class utype> bool
    isidentity(const xmatrix<utype>& a) {
      //ME20200725 - changed exit to return false
      if(a.rows!=a.cols) return false;
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          if(i-a.lrows+1!=j-a.lcols+1){  // i != j
            if(aurostd::abs(a[i][j]) >_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
          }else{
            if(aurostd::abs(1.0-a[i][j]) >_AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
          }
      return TRUE;
    }
  //CO START

  // namespace aurostd
  template<class utype> bool
    isdiagonal(const xmatrix<utype>& a,const utype& _eps_) { //DX20171025
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          if(i-a.lrows+1!=j-a.lcols+1)  // i != j
            if(aurostd::abs(a[i][j]) >_eps_) return FALSE;
      return TRUE;
    }

  // namespace aurostd
  template<class utype> bool
    issymmetric(const xmatrix<utype>& a) {
      //ME20200725 - changed exit to return false
      if(a.rows!=a.cols) return false;
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          if(aurostd::abs(a[i][j]-a[j][i])> _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
      return TRUE;
    }

  // namespace aurostd
  template<class utype> bool
    isantisymmetric(const xmatrix<utype>& a) {
      //ME20200725 - changed exit to return false
      if(a.rows!=a.cols) return false;
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          if(aurostd::abs(a[i][j]-(-a[j][i]))> _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
      return TRUE;
    }

  // namespace aurostd
  template<class utype> bool
    ishermitian(const xmatrix<utype>& a) {
      //ME20200725 - changed exit to return false
      if(a.rows!=a.cols) return false;
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++) {
          if(aurostd::abs(a[i][j]-aurostd::conj(a[j][i])) > _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
        }
      return TRUE;
    }

  // namespace aurostd
  template<class utype> bool
    isantihermitian(const xmatrix<utype>& a) {
      //ME20200725 - changed exit to return false
      if(a.rows!=a.cols) return false;
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++) {
          if(aurostd::abs(a[i][j]-(-aurostd::conj(a[j][i]))) > _AUROSTD_XMATRIX_TOLERANCE_IDENTITY_) return FALSE;
        }
      return TRUE;
    }

}

// ****************************************************************************
// ------------------------------------------------------ xmatrix construction
namespace aurostd {
  /// @brief reshape xvector to xmatrix with given dimensions
  /// @param vec input vector
  /// @param rows row count
  /// @param cols column count
  /// @return created matrix
  /// @authors
  /// @mod{ME,20210921,created function}
  /// @mod{HE,20220915,migrated to direct reshape}
  /// @note When rows or cols is set to 0, it is automatically calculated.
  /// This is also true when cols is not given.
  template<class utype>
  xmatrix<utype> reshape(const xvector<utype> &vec, int rows, int cols) {
    xmatrix<utype> mat = reshape(vec);
    mat.reshape(rows, cols);
    return mat;
  }

  /// @brief reshape xvector to xmatrix column
  /// @param vec input vector
  /// @return created xmatrix
  /// @authors
  /// @mod{ME,20190617,created function}
  /// @mod{HE,20220915,using general stacking function}
  template<class utype>
  xmatrix<utype> reshape(const xvector<utype> &vec) {
    xmatrix<utype> mat(vec.rows, 1);
    for (int i = 1; i <= vec.rows; i++)
      mat(i-vec.lrows+1, 1) = vec(i);
    return mat;
  }

  /// @brief reshape xvector to xmatrix row
  /// @param vec input vector
  /// @return created xmatrix
  /// @authors
  /// @mod{ME,20190617,created function}
  /// @mod{HE,20220915,using general stacking function}
  template<class utype>
  xmatrix<utype> reshape_rows(const xvector<utype> &vec) {
    xmatrix<utype> mat = reshape(vec);
    mat.reshape(1);
    return mat;
  }

  /// @brief stack xvectors vertically into a xmatrix
  /// @param vec_list list of xvectors
  /// @return generated xmatrix
  /// @authors
  /// @mod{HE,20220915,created function}
  /// @code{.cpp}
  /// xvector<double> xv = {1.5, 3.0, 4.5};
  /// xmatrix<double> smat = aurostd::vstack(vector<xvector<double>>({xv, xv, xv}));
  /// cout << smat << endl;
  /// // 1.5 1.5 1.5
  /// // 3.0 3.0 3.0
  /// // 4.5 4.5 4.5
  /// @endcode
  template<class utype>
  xmatrix<utype> vstack(const vector <xvector<utype>> &vec_list) {
    if (vec_list.empty()) return xmatrix<utype>();
    int rows = vec_list[0].rows;
    int cols = (int) vec_list.size();
    for (xvector<utype> const &vec: vec_list) {
      if (rows != vec.rows) {
        stringstream message;
        message << "vectors must have the same dimensions " << rows << " != " << vec.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
    }
    xmatrix<utype> mat(rows, cols);
    for (int i = 1; i <= rows; i++) {
      for (int j = 1; j <= cols; j++) {
        mat[i][j] = vec_list[j-1][i];
      }
    }
    return mat;
  }



  /// @brief Stack xvectors horizontally into a xmatrix
  /// @param vec_list list of xvectors
  /// @return generated xmatrix
  /// @authors
  /// @mod{HE,20220915,created function}
  /// @note use list initializing to call this function
  /// @code{.cpp}
  /// xvector<double> xv = {1.5, 3.0, 4.5};
  /// xmatrix<double> smat = aurostd::hstack(vector<xvector<double>>({xv, xv, xv}));
  /// cout << smat << endl;
  /// // 1.5 3.0 4.5
  /// // 1.5 3.0 4.5
  /// // 1.5 3.0 4.5
  /// @endcode
  template<class utype>
  xmatrix<utype> hstack(const vector <xvector<utype>> &vec_list) {
    if (vec_list.empty()) return xmatrix<utype>();
    int cols = vec_list[0].rows;
    int rows = (int) vec_list.size();
    for (xvector<utype> const &vec: vec_list) {
      if (cols != vec.rows) {
        stringstream message;
        message << "vectors must have the same dimensions " << rows << " != " << vec.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
    }
    xmatrix<utype> mat(rows, cols);
    for (int i = 1; i <= rows; i++) {
      for (int j = 1; j <= cols; j++) {
        mat[i][j] = vec_list[i-1][j];
      }
    }
    return mat;
  }

}

// ****************************************************************************
// -------------------------------------------------------------- xmatrix example types

namespace aurostd {  // namespace aurostd
  //20171008 - CO
  //doesn't have to be square like identity
  template<class utype> xmatrix<utype>
    eye(int nrh,int nch,int nrl,int ncl) __xprototype { //CO20190520
      if(nch==AUROSTD_MAX_INT){nch=nrh;}  //eye(3)==eye(3,3)
      xmatrix<utype> a(nrh,nch,nrl,ncl);
      //[CO20200106 - doesn't work for lrows!=lcols]for (int i=a.lrows;i<=a.urows && i<=a.ucols;i++){a[i][i] = (utype)1;} //ME20200106
      int i=0,j=0;
      for (i=a.lrows;i<=a.urows;i++){
        for (j=a.lcols;j<=a.ucols;j++){
          if(i==j){a[i][j]=(utype)1;}
        }
      }
      return a;
    }
}

namespace aurostd {  // namespace aurostd
  //20171008 - CO
  template<class utype> xmatrix<utype>
    ones_xm(int nrh,int nch,int nrl,int ncl) __xprototype { //CO20190520
      xmatrix<utype> a(nrh,nch,nrl,ncl);
      for (int i=a.lrows;i<=a.urows;i++){
        for (int j=a.lcols;j<=a.ucols;j++){
          a[i][j]=(utype)1;
        }
      }
      return a;
    }
}

// ****************************************************************************
// -------------------------------------------------------------- xmatrix casts

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to long double
    xmatrix<long double> xlongdouble(const xmatrix<utype> &a) __xprototype {
      xmatrix<long double> c(a.urows,a.ucols,a.lrows,a.lcols);
      for (int i=a.lrows;i<=a.urows;i++)
        for (int j=a.lcols;j<=a.ucols;j++)
          c[i][j]=(long double) a[i][j];
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to double
    xmatrix<double> xdouble(const xmatrix<utype> &a) __xprototype {
      xmatrix<double> c(a.urows,a.ucols,a.lrows,a.lcols);
      for (int i=a.lrows;i<=a.urows;i++)
        for (int j=a.lcols;j<=a.ucols;j++)
          c[i][j]=(double) a[i][j];
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to float
    xmatrix<float> xfloat(const xmatrix<utype> &a) __xprototype {
      xmatrix<float> c(a.urows,a.ucols,a.lrows,a.lcols);
      for (int i=a.lrows;i<=a.urows;i++)
        for (int j=a.lcols;j<=a.ucols;j++)
          c[i][j]=(float) a[i][j];
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to long int
    xmatrix<long int> xlongint(const xmatrix<utype> &a) __xprototype {
      xmatrix<long int> c(a.urows,a.ucols,a.lrows,a.lcols);
      for (int i=a.lrows;i<=a.urows;i++)
        for (int j=a.lcols;j<=a.ucols;j++)
          c[i][j]=(long int) a[i][j];
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to int
    xmatrix<int> xint(const xmatrix<utype> &a) __xprototype {
      xmatrix<int> c(a.urows,a.ucols,a.lrows,a.lcols);
      for (int i=a.lrows;i<=a.urows;i++)
        for (int j=a.lcols;j<=a.ucols;j++)
          c[i][j]=(int) a[i][j];
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to char
    xmatrix<char> xchar(const xmatrix<utype> &a) __xprototype {
      xmatrix<char> c(a.urows,a.ucols,a.lrows,a.lcols);
      for (int i=a.lrows;i<=a.urows;i++)
        for (int j=a.lcols;j<=a.ucols;j++)
          c[i][j]=(char) a[i][j];
      return c;
    }
}

namespace aurostd {                   // conversion to vector<vector<utype> >
  template<class utype> vector<vector<utype> >
    xmatrix2vectorvector(const xmatrix<utype>& xmat) __xprototype {
      int isize=xmat.rows,jsize=xmat.cols;
      // vector<vector<utype> > mat; vector<utype> v; mat=vector<vector<utype> > (m,v); // by hand
      // vector<vector<utype> > vectorvector(isize,jsize);              // WORKS WITH gcc/g++ 4.2 and 4.1
      vector<vector<utype> > vectorvector(isize,vector<utype>(jsize));  // WORKS WITH gcc/g++ 4.3

      for(int i=0;i<isize;i++)
        for(int j=0;j<jsize;j++)
          vectorvector[i][j]=xmat(i+xmat.lrows,j+xmat.lcols);
      return vectorvector;
    }
}

//SD20220504 - Rewritten vectorvector2xmatrix to be more like vector2xvector
namespace aurostd {                   // conversion to xmatrix<utype>
  template<class utype> xmatrix<utype>
    vectorvector2xmatrix(const vector<vector<utype> >& mat, int lrows, int lcols) __xprototype {
      int isize=mat.size(),jsize=mat[0].size();
      xmatrix<utype> xmat(isize+lrows-1,jsize+lcols-1);
      for(int i=lrows;i<=isize+lrows-1;i++)
        for(int j=lcols;j<=jsize+lcols-1;j++)
          xmat[i][j]=mat[i-lrows][j-lcols];
      return xmat;
    }
  template<class utype> xmatrix<utype>
    vectorvector2xmatrix(const vector<vector<string> >& mat, int lrows, int lcols) __xprototype {
      int isize=mat.size(),jsize=mat[0].size();
      xmatrix<utype> xmat(isize+lrows-1,jsize+lcols-1);
      for(int i=lrows;i<=isize+lrows-1;i++)
        for(int j=lcols;j<=jsize+lcols-1;j++)
          xmat[i][j]=aurostd::string2utype<utype>(mat[i-lrows][j-lcols]);
      return xmat;
    }
}

//CO20191201
namespace aurostd {                   // conversion from xmatrix<int> to xmatrix<double>
  template<class utype> xmatrix<double>
    xmatrixutype2double(const xmatrix<utype>& a){ //CO20191201
      xmatrix<double> b(a.urows,a.ucols,a.lrows,a.lcols);
      int i=0,j=0;
      for(i=a.lrows;i<=a.urows;i++){
        for(j=a.lcols;j<=a.ucols;j++){
          b[i][j]=(double)a[i][j];
        }
      }
      return b;
    }
}

//CO20191201
namespace aurostd {                   // conversion from xmatrix<int> to xmatrix<double>
  template<class utype> xmatrix<utype>
    xmatrixdouble2utype(const xmatrix<double>& a,bool check_int){  //CO20191201
      xmatrix<utype> b(a.urows,a.ucols,a.lrows,a.lcols);
      int i=0,j=0;
      if(check_int){
        for(i=a.lrows;i<=a.urows;i++){
          for(j=a.lcols;j<=a.ucols;j++){
            if(!isinteger(a[i][j])){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::xmatrixdouble2utype():","non-integer found ["+aurostd::utype2string(a[i][j])+"]",_INPUT_ILLEGAL_);}
          }
        }
      }
      for(i=a.lrows;i<=a.urows;i++){
        for(j=a.lcols;j<=a.ucols;j++){
          b[i][j]=(utype)nint(a[i][j]);  //nint is for safety
        }
      }
      return b;
    }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                               // function reset xmatrix<>
    void xmatrix<utype>::reset(void) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      cout<<"M -> function reset: "
        <<" lrows="<<lrows<<" urows="<<urows<<" lcols="<<lcols<<" ucols="<<ucols<<endl;
#endif
      for(int i=lrows;i<=urows;i++)
        for(int j=lcols;j<=ucols;j++)
          corpus[i][j]=(utype) 0.0;
    }
  template<class utype>                               // function reset xmatrix<>
    void reset(xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      cout<<"M -> function reset: "
        <<" a.lrows="<<a.lrows<<" a.urows="<<a.urows<<" a.lcols="<<a.lcols<<" a.ucols="<<a.ucols<<endl;
#endif
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          a[i][j]=(utype) 0.0;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                               // function clear xmatrix<>
    void xmatrix<utype>::clear(void) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      cout<<"M -> function clear: "
        <<" lrows="<<lrows<<" urows="<<urows<<" lcols="<<lcols<<" ucols="<<ucols<<endl;
#endif
      reset();
    }
  template<class utype>                               // function clear xmatrix<>
    void clear(xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      cout<<"M -> function clear: "
        <<" a.lrows="<<a.lrows<<" a.urows="<<a.urows<<" a.lcols="<<a.lcols<<" a.ucols="<<a.ucols<<endl;
#endif
      reset(a);
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function set xmatrix<>
    void xmatrix<utype>::set(const utype& s) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      cout<<"M -> function set: "
        <<" lrows="<<lrows<<" urows="<<urows<<" lcols="<<lcols<<" ucols="<<ucols<<endl;
#endif
      for(int i=lrows;i<=urows;i++)
        for(int j=lcols;j<=ucols;j++)
          corpus[i][j]=(utype) s;
    }
  template<class utype>                                 // function set xmatrix<>
    void set(xmatrix<utype>& a,const utype& s) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      cout<<"M -> function set: "
        <<" a.lrows="<<a.lrows<<" a.urows="<<a.urows<<" a.lcols="<<a.lcols<<" a.ucols="<<a.ucols<<endl;
#endif
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          a[i][j]=(utype) s;
    }
}

// ----------------------------------------------------------------------------
//[SD20220616 - OBSOLETE]namespace aurostd {  // namespace aurostd
//[SD20220616 - OBSOLETE]  template<class utype>                           // function vector<xmatrix<>>
//[SD20220616 - OBSOLETE]    xvector<utype> vector(const xmatrix<utype>& a) {
//[SD20220616 - OBSOLETE]      int n=(a.rows*a.cols);
//[SD20220616 - OBSOLETE]      xvector<utype> c(1,n);
//[SD20220616 - OBSOLETE]      for(int i=0;i<n;i++) {
//[SD20220616 - OBSOLETE]        c[i+1]=(utype) a(int(i/a.cols)+a.lrows,mod(i,a.cols)+a.lcols);
//[SD20220616 - OBSOLETE]      }
//[SD20220616 - OBSOLETE]      return c;
//[SD20220616 - OBSOLETE]    }
//[SD20220616 - OBSOLETE]}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function det xmatrix<>
    utype det(const xmatrix<utype>& a) {
      /* returns the determinant **/
      if(!a.issquare){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::det()","a must be square",_INPUT_ILLEGAL_);}
      if(a.lrows!=1 || a.lcols!=1){xmatrix<utype> b=a;b.shift(1,1);return det(b);}
      int size=a.rows;
      //  cerr << "DET CALL size="<<size<< endl;
      if(size==1) { return (utype) a[1][1]; }
      if(size==2) { return (utype) a[1][1]*a[2][2]-a[1][2]*a[2][1]; }
      if(size==3) {
        return (utype) (a[1][1]*a[2][2]*a[3][3]+
            a[1][2]*a[2][3]*a[3][1]+
            a[1][3]*a[2][1]*a[3][2]-
            a[1][3]*a[2][2]*a[3][1]-
            a[1][2]*a[2][1]*a[3][3]-
            a[1][1]*a[2][3]*a[3][2]); }
      if(size==4) {
        return (utype) (a[1][4]*a[2][3]*a[3][2]*a[4][1]-a[1][3]*a[2][4]*a[3][2]*a[4][1]-a[1][4]*a[2][2]*a[3][3]*a[4][1]+a[1][2]*a[2][4]*a[3][3]*a[4][1]+
            a[1][3]*a[2][2]*a[3][4]*a[4][1]-a[1][2]*a[2][3]*a[3][4]*a[4][1]-a[1][4]*a[2][3]*a[3][1]*a[4][2]+a[1][3]*a[2][4]*a[3][1]*a[4][2]+
            a[1][4]*a[2][1]*a[3][3]*a[4][2]-a[1][1]*a[2][4]*a[3][3]*a[4][2]-a[1][3]*a[2][1]*a[3][4]*a[4][2]+a[1][1]*a[2][3]*a[3][4]*a[4][2]+
            a[1][4]*a[2][2]*a[3][1]*a[4][3]-a[1][2]*a[2][4]*a[3][1]*a[4][3]-a[1][4]*a[2][1]*a[3][2]*a[4][3]+a[1][1]*a[2][4]*a[3][2]*a[4][3]+
            a[1][2]*a[2][1]*a[3][4]*a[4][3]-a[1][1]*a[2][2]*a[3][4]*a[4][3]-a[1][3]*a[2][2]*a[3][1]*a[4][4]+a[1][2]*a[2][3]*a[3][1]*a[4][4]+
            a[1][3]*a[2][1]*a[3][2]*a[4][4]-a[1][1]*a[2][3]*a[3][2]*a[4][4]-a[1][2]*a[2][1]*a[3][3]*a[4][4]+a[1][1]*a[2][2]*a[3][3]*a[4][4]);  }
      utype out=(utype) 0;
      if(size>=5) {
        xmatrix<utype> b(size-1,size-1);
        for(int j=1;j<=size;j++) {
          for(int ib=1;ib<=size-1;ib++)                                         // make sub
            for(int jb=1;jb<=size-1;jb++)                                       // make sub
              if(jb<j) { b[ib][jb]=a[ib+1][jb]; } else { b[ib][jb]=a[ib+1][jb+1]; }  // make sub --- FASTER
          if(_isodd(j)) { out+=a[1][j]*det(b); } else { out-=a[1][j]*det(b); }         // get part --- FASTER
          //if(jb<j)  { b(ib,jb)=a(ib+1,jb); } else {  b(ib,jb)=a(ib+1,jb+1); }        // make sub --- SLOWER
          //if(_isodd(j)) { out+=a(1,j)*det(b); } else { out-=a(1,j)*det(b); }         // get part --- SLOWER
        }
        return (utype) out;
      }
      return (utype) out;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                       // function determinant xmatrix<>
    utype determinant(const xmatrix<utype>& a) {
      return (utype) det(a);
    }
}

// ----------------------------------------------------------------------------
// Cayley-Menger Determinant
// http://mathworld.wolfram.com/Cayley-MengerDeterminant.html
//CO20180515
namespace aurostd { // namespace aurostd
  template<class utype> utype
    CMdet(const xmatrix<utype>& B){ //Cayley-Menger Determinant for simplex content
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      if(!B.issquare){
        string message = "Only defined for square xmatrices";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      xmatrix<utype> B_hat=ones_xm<utype>(B.urows+1,B.ucols+1,B.lrows,B.lcols); //CO20190520
      B_hat(B.lrows,B.lcols)=(utype)0;
      for(int row=B.lrows;row<=B.urows;row++){
        for(int col=B.lcols;col<=B.ucols;col++){
          B_hat(row+1,col+1)=B(row,col);
        }
      }
      utype detB_hat=det(B_hat);
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " B_hat=" << endl; cerr << B_hat << endl;
        cerr << __AFLOW_FUNC__ << " det(B_hat)=" << detB_hat << endl;
      }
      return detB_hat;
    }
}

// ----------------------------------------------------------------------------
// getRotationMatrix3D() rotating a onto b
// https://math.stackexchange.com/questions/20180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
//CO20190324
namespace aurostd { // namespace aurostd
  template<class utype> xmatrix<utype>
    getRotationMatrix3D(const xvector<utype>& a,const xvector<utype>& b){
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      //tests of stupidity
      if(a.rows!=3){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"a.rows!=3",_INPUT_ILLEGAL_);}
      if(b.rows!=3){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"b.rows!=3",_INPUT_ILLEGAL_);}
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " a=" << a << endl;
        cerr << __AFLOW_FUNC__ << " b=" << b << endl;
      }
      xmatrix<utype> R=aurostd::eye<utype>(3,3); //CO20190520

      //trivial cases
      if(aurostd::isequal(a,b,(utype)_ZERO_TOL_)){return R;}
      if(aurostd::isequal(a,-b,(utype)_ZERO_TOL_)){R(1,1)=R(2,2)=R(3,3)=-1;return R;}

      //non-trivial case
      xvector<utype> v=aurostd::vector_product(a,b); //requires first index 1
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " v=" << v << endl;}
      utype s=aurostd::modulus(v);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " s=" << s << endl;}
      utype c=aurostd::scalar_product(a,b);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " c=" << c << endl;}
      xmatrix<utype> vx(3,3);
      vx[2][1]=v[3];
      vx[3][1]=-v[2];
      vx[1][2]=-v[3];
      vx[3][2]=v[1];
      vx[1][3]=v[2];
      vx[2][3]=-v[1];
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vx=" << endl;cerr << vx << endl;}
      if(!aurostd::isequal(s,(utype)0,(utype)_ZERO_TOL_)){R+=vx+vx*vx*((utype)1-c)/(s*s);}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " R=" << endl;cerr << R << endl;}
      return R;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> void                       // minor submatrix - IN PLACE (faster for big routines)
    submatrixInPlace(const xmatrix<utype>& a,xmatrix<utype>& b,int irow,int jcol) {
      if(irow>=a.lrows && irow<=a.urows && jcol>=a.lcols && jcol<=a.ucols) {
        //[CO20191201 - OBSOLETE]xmatrix<utype> b(a.urows-1,a.ucols-1,a.lrows,a.lcols);
        if(b.lrows!=a.lrows || b.lcols!=a.lcols || b.urows!=a.urows-1 || b.ucols!=b.ucols-1){ //CO20191201
          xmatrix<utype> c(a.urows-1,a.ucols-1,a.lrows,a.lcols);
          b=c;
        }
        for(int i=a.lrows;i<=a.urows;i++)
          for(int j=a.lcols;j<=a.ucols;j++) {
            if(i<irow && j<jcol) b[i][j]=a[i][j];
            if(i>irow && j<jcol) b[i-1][j]=a[i][j];
            if(i<irow && j>jcol) b[i][j-1]=a[i][j];
            if(i>irow && j>jcol) b[i-1][j-1]=a[i][j];
          }
        //[CO20191201 - OBSOLETE]return b;
        return;
      }
      if((irow<a.lrows || irow>a.urows) && jcol>=a.lcols && jcol<=a.ucols) {
        //[CO20191201 - OBSOLETE]xmatrix<utype> b(a.urows,a.ucols-1,a.lrows,a.lcols);
        if(b.lrows!=a.lrows || b.lcols!=a.lcols || b.urows!=a.urows || b.ucols!=b.ucols-1){ //CO20191201
          xmatrix<utype> c(a.urows,a.ucols-1,a.lrows,a.lcols);
          b=c;
        }
        for(int i=a.lrows;i<=a.urows;i++)
          for(int j=a.lcols;j<=a.ucols;j++) {
            if(j<jcol) b[i][j]=a[i][j];
            if(j>jcol) b[i][j-1]=a[i][j];
          }
        //[CO20191201 - OBSOLETE]return b;
        return;
      }
      if(irow>=a.lrows && irow<=a.urows && (jcol<a.lcols || jcol>a.ucols)) {
        //[CO20191201 - OBSOLETE]xmatrix<utype> b(a.urows-1,a.ucols,a.lrows,a.lcols);
        if(b.lrows!=a.lrows || b.lcols!=a.lcols || b.urows!=a.urows-1 || b.ucols!=b.ucols){ //CO20191201
          xmatrix<utype> c(a.urows-1,a.ucols,a.lrows,a.lcols);
          b=c;
        }
        for(int i=a.lrows;i<=a.urows;i++)
          for(int j=a.lcols;j<=a.ucols;j++) {
            if(i<irow) b[i][j]=a[i][j];
            if(i>irow) b[i-1][j]=a[i][j];
          }
        //[CO20191201 - OBSOLETE]return b;
        return;
      }
      //[CO20191201 - OBSOLETE]return a;
      b=a;  //CO20191201
    }
  template<class utype> xmatrix<utype>                       // minor submatrix - IN PLACE
    submatrix(const xmatrix<utype>& a,int irow,int jcol) {
      if(irow>=a.lrows && irow<=a.urows && jcol>=a.lcols && jcol<=a.ucols) {
        xmatrix<utype> b(a.urows-1,a.ucols-1,a.lrows,a.lcols);
        submatrixInPlace(a,b,irow,jcol);
        return b;
      }
      if((irow<a.lrows || irow>a.urows) && jcol>=a.lcols && jcol<=a.ucols) {
        xmatrix<utype> b(a.urows,a.ucols-1,a.lrows,a.lcols);
        submatrixInPlace(a,b,irow,jcol);
        return b;
      }
      if(irow>=a.lrows && irow<=a.urows && (jcol<a.lcols || jcol>a.ucols)) {
        xmatrix<utype> b(a.urows-1,a.ucols,a.lrows,a.lcols);
        submatrixInPlace(a,b,irow,jcol);
        return b;
      }
      return a;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> utype                       // minor submatrix
    minordet(const xmatrix<utype>& a,const int& irow,const int& jcol) {
      return det(submatrix(a,irow,jcol));
    }
  template<class utype> utype                       // minor submatrix
    minordeterminant(const xmatrix<utype>& a,const int& irow,const int& jcol) {
      return determinant(submatrix(a,irow,jcol));
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function inverse xmatrix<>
    void adjointInPlace(const xmatrix<utype>& a,xmatrix<utype>& b) { //CO20191201
      //inspired by https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
      if(!a.issquare){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::inverseByAdjoint()","a must be square",_INPUT_ILLEGAL_);}
      b=a;
      xmatrix<utype> submat(a.urows-1,a.ucols-1,a.lrows,a.lcols);
      int i=0,j=0;
      for(i=a.lrows;i<=a.urows;i++){
        for(j=a.lcols;j<=a.ucols;j++){
          submatrixInPlace(a,submat,i,j);
          b[i][j]=(utype)aurostd::powint(-1,i+j)*det(submat);
        }
      }
      traspInPlace(b);
    }

  template<class utype>                                 // function inverse xmatrix<>
    xmatrix<utype> adjoint(const xmatrix<utype>& a) { //CO20191201
      xmatrix<utype> b;
      adjointInPlace(a,b);
      return b;
    }

  template<class utype>                                 // function inverse xmatrix<>
    xmatrix<utype> inverseByAdjoint(const xmatrix<utype>& a) {return (utype)1.0/det(a) * adjoint(a);} //CO20191201

  /// @brief calculates inverse of a square matrix
  ///
  /// @param a square matrix
  ///
  /// @return inverse of the square matrix
  ///
  /// @authors
  /// @mod{CO,20191201,created function}
  /// @mod{SD,20220427,added inverse by QR decomposition}
  ///
  /// @note The function will first try to find the inverse using the adjoint and if that fails
  /// then using QR decomposition
  template<class utype>                                 // function inverse xmatrix<>
    xmatrix<utype> inverse(const xmatrix<utype>& a) {
      // returns the inverse
      if(!a.issquare){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"a must be square",_INPUT_ILLEGAL_);}
      if(a.lrows!=1 || a.lcols!=1){
        xmatrix<utype> b(a);
        b.shift(1,1);
        b=inverse(b);
        b.shift(a.lrows,a.lcols);
        return b;
      }
      int size=a.rows;
      xmatrix<utype> b(a.rows,a.cols);
      //  cerr << "DET CALL size="<<size<< endl;
      utype adet=det(a);
      if(adet==(utype) 0)  {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"singular matrix",_INPUT_ILLEGAL_);}
      if(size==1) {b[1][1]=(utype)1/a[1][1]; return b;}
      if(size==2) { //CO20191201
        b[1][1]=a[2][2]/adet;b[1][2]=-a[1][2]/adet;
        b[2][1]=-a[2][1]/adet;b[2][2]=a[1][1]/adet;
        return b;
      }
      if(size==3) {
        //[CO20191201 - already calculated above]adet=det(a);
        //b[1][1]=(+a[2][2]*a[3][3]-a[2][3]*a[3][2])/adet;   // with fast index []
        //b[1][2]=(-a[1][2]*a[3][3]+a[1][3]*a[3][2])/adet;   // with fast index []
        //b[1][3]=(+a[1][2]*a[2][3]-a[1][3]*a[2][2])/adet;   // with fast index []
        //b[2][1]=(-a[2][1]*a[3][3]+a[2][3]*a[3][1])/adet;   // with fast index []
        //b[2][2]=(+a[1][1]*a[3][3]-a[1][3]*a[3][1])/adet;   // with fast index []
        //b[2][3]=(-a[1][1]*a[2][3]+a[1][3]*a[2][1])/adet;   // with fast index []
        //b[3][1]=(+a[2][1]*a[3][2]-a[2][2]*a[3][1])/adet;   // with fast index []
        //b[3][2]=(-a[1][1]*a[3][2]+a[1][2]*a[3][1])/adet;   // with fast index []
        //b[3][3]=(+a[1][1]*a[2][2]-a[1][2]*a[2][1])/adet;   // with fast index []
        b(1,1)=(+a(2,2)*a(3,3)-a(2,3)*a(3,2))/adet;   // with slow index ()
        b(1,2)=(-a(1,2)*a(3,3)+a(1,3)*a(3,2))/adet;   // with slow index ()
        b(1,3)=(+a(1,2)*a(2,3)-a(1,3)*a(2,2))/adet;   // with slow index ()
        b(2,1)=(-a(2,1)*a(3,3)+a(2,3)*a(3,1))/adet;   // with slow index ()
        b(2,2)=(+a(1,1)*a(3,3)-a(1,3)*a(3,1))/adet;   // with slow index ()
        b(2,3)=(-a(1,1)*a(2,3)+a(1,3)*a(2,1))/adet;   // with slow index ()
        b(3,1)=(+a(2,1)*a(3,2)-a(2,2)*a(3,1))/adet;   // with slow index ()
        b(3,2)=(-a(1,1)*a(3,2)+a(1,2)*a(3,1))/adet;   // with slow index ()
        b(3,3)=(+a(1,1)*a(2,2)-a(1,2)*a(2,1))/adet;   // with slow index ()
        return b;
      }
      if(size==4) {
        //[CO20191201 - already calculated above]adet=det(a);
        //b[1][1]=(+a[2][2]*(a[3][3]*a[4][4]-a[3][4]*a[4][3])-a[2][3]*(a[3][2]*a[4][4]-a[3][4]*a[4][2])+a[2][4]*(a[3][2]*a[4][3]-a[3][3]*a[4][2]))/adet;   // with fast index []
        //b[2][1]=(-a[2][1]*(a[3][3]*a[4][4]-a[3][4]*a[4][3])+a[2][3]*(a[3][1]*a[4][4]-a[3][4]*a[4][1])-a[2][4]*(a[3][1]*a[4][3]-a[3][3]*a[4][1]))/adet;   // with fast index []
        //b[3][1]=(+a[2][1]*(a[3][2]*a[4][4]-a[3][4]*a[4][2])-a[2][2]*(a[3][1]*a[4][4]-a[3][4]*a[4][1])+a[2][4]*(a[3][1]*a[4][2]-a[3][2]*a[4][1]))/adet;   // with fast index []
        //b[4][1]=(-a[2][1]*(a[3][2]*a[4][3]-a[3][3]*a[4][2])+a[2][2]*(a[3][1]*a[4][3]-a[3][3]*a[4][1])-a[2][3]*(a[3][1]*a[4][2]-a[3][2]*a[4][1]))/adet;   // with fast index []
        //b[1][2]=(-a[1][2]*(a[3][3]*a[4][4]-a[3][4]*a[4][3])+a[1][3]*(a[3][2]*a[4][4]-a[3][4]*a[4][2])-a[1][4]*(a[3][2]*a[4][3]-a[3][3]*a[4][2]))/adet;   // with fast index []
        //b[2][2]=(+a[1][1]*(a[3][3]*a[4][4]-a[3][4]*a[4][3])-a[1][3]*(a[3][1]*a[4][4]-a[3][4]*a[4][1])+a[1][4]*(a[3][1]*a[4][3]-a[3][3]*a[4][1]))/adet;   // with fast index []
        //b[3][2]=(-a[1][1]*(a[3][2]*a[4][4]-a[3][4]*a[4][2])+a[1][2]*(a[3][1]*a[4][4]-a[3][4]*a[4][1])-a[1][4]*(a[3][1]*a[4][2]-a[3][2]*a[4][1]))/adet;   // with fast index []
        //b[4][2]=(+a[1][1]*(a[3][2]*a[4][3]-a[3][3]*a[4][2])-a[1][2]*(a[3][1]*a[4][3]-a[3][3]*a[4][1])+a[1][3]*(a[3][1]*a[4][2]-a[3][2]*a[4][1]))/adet;   // with fast index []
        //b[1][3]=(+a[1][2]*(a[2][3]*a[4][4]-a[2][4]*a[4][3])-a[1][3]*(a[2][2]*a[4][4]-a[2][4]*a[4][2])+a[1][4]*(a[2][2]*a[4][3]-a[2][3]*a[4][2]))/adet;   // with fast index []
        //b[2][3]=(-a[1][1]*(a[2][3]*a[4][4]-a[2][4]*a[4][3])+a[1][3]*(a[2][1]*a[4][4]-a[2][4]*a[4][1])-a[1][4]*(a[2][1]*a[4][3]-a[2][3]*a[4][1]))/adet;   // with fast index []
        //b[3][3]=(+a[1][1]*(a[2][2]*a[4][4]-a[2][4]*a[4][2])-a[1][2]*(a[2][1]*a[4][4]-a[2][4]*a[4][1])+a[1][4]*(a[2][1]*a[4][2]-a[2][2]*a[4][1]))/adet;   // with fast index []
        //b[4][3]=(-a[1][1]*(a[2][2]*a[4][3]-a[2][3]*a[4][2])+a[1][2]*(a[2][1]*a[4][3]-a[2][3]*a[4][1])-a[1][3]*(a[2][1]*a[4][2]-a[2][2]*a[4][1]))/adet;   // with fast index []
        //b[1][4]=(-a[1][2]*(a[2][3]*a[3][4]-a[2][4]*a[3][3])+a[1][3]*(a[2][2]*a[3][4]-a[2][4]*a[3][2])-a[1][4]*(a[2][2]*a[3][3]-a[2][3]*a[3][2]))/adet;   // with fast index []
        //b[2][4]=(+a[1][1]*(a[2][3]*a[3][4]-a[2][4]*a[3][3])-a[1][3]*(a[2][1]*a[3][4]-a[2][4]*a[3][1])+a[1][4]*(a[2][1]*a[3][3]-a[2][3]*a[3][1]))/adet;   // with fast index []
        //b[3][4]=(-a[1][1]*(a[2][2]*a[3][4]-a[2][4]*a[3][2])+a[1][2]*(a[2][1]*a[3][4]-a[2][4]*a[3][1])-a[1][4]*(a[2][1]*a[3][2]-a[2][2]*a[3][1]))/adet;   // with fast index []
        //b[4][4]=(+a[1][1]*(a[2][2]*a[3][3]-a[2][3]*a[3][2])-a[1][2]*(a[2][1]*a[3][3]-a[2][3]*a[3][1])+a[1][3]*(a[2][1]*a[3][2]-a[2][2]*a[3][1]))/adet;   // with fast index []
        b(1,1)=(+a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))/adet;   // with slow index ()
        b(2,1)=(-a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))/adet;   // with slow index ()
        b(3,1)=(+a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))/adet;   // with slow index ()
        b(4,1)=(-a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))-a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))/adet;   // with slow index ()
        b(1,2)=(-a(1,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(1,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(1,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))/adet;   // with slow index ()
        b(2,2)=(+a(1,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(1,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(1,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))/adet;   // with slow index ()
        b(3,2)=(-a(1,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(1,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(1,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))/adet;   // with slow index ()
        b(4,2)=(+a(1,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(1,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(1,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))/adet;   // with slow index ()
        b(1,3)=(+a(1,2)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))-a(1,3)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))+a(1,4)*(a(2,2)*a(4,3)-a(2,3)*a(4,2)))/adet;   // with slow index ()
        b(2,3)=(-a(1,1)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))+a(1,3)*(a(2,1)*a(4,4)-a(2,4)*a(4,1))-a(1,4)*(a(2,1)*a(4,3)-a(2,3)*a(4,1)))/adet;   // with slow index ()
        b(3,3)=(+a(1,1)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))-a(1,2)*(a(2,1)*a(4,4)-a(2,4)*a(4,1))+a(1,4)*(a(2,1)*a(4,2)-a(2,2)*a(4,1)))/adet;   // with slow index ()
        b(4,3)=(-a(1,1)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))+a(1,2)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))-a(1,3)*(a(2,1)*a(4,2)-a(2,2)*a(4,1)))/adet;   // with slow index ()
        b(1,4)=(-a(1,2)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))+a(1,3)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))-a(1,4)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)))/adet;   // with slow index ()
        b(2,4)=(+a(1,1)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))-a(1,3)*(a(2,1)*a(3,4)-a(2,4)*a(3,1))+a(1,4)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)))/adet;   // with slow index ()
        b(3,4)=(-a(1,1)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))+a(1,2)*(a(2,1)*a(3,4)-a(2,4)*a(3,1))-a(1,4)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))/adet;   // with slow index ()
        b(4,4)=(+a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))/adet;   // with slow index ()
        return b;
      }
      //CO20191201 - GaussJordan() is INCREDIBLY unstable (division by small numbers over and over again)
      //use inverseByAdjoint, which waits until the end to divide by a (hopefully) bigger number
      //[CO20191201 - OBSOLETE]// if everything fails move to GaussJordan
      //[CO20191201 - OBSOLETE]b=a;
      //[CO20191201 - OBSOLETE]xmatrix<utype> B(a.rows,a.cols);
      //[CO20191201 - OBSOLETE]GaussJordan(b,B);
      //    if(size>=6) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xmatrix<utype>::inverse: " << size << "x" << size << " not written yet" << endl;}
      //[CO20191201 - OBSOLETE]return b;
      //SD20220427 - start with finding the inverse by adjoint since an answer is guaranteed, if the answer is wrong, then go on to other methods
      xmatrix<utype> id = aurostd::identity(a);
      try { //SD20220427 - we need a try block because LUP and QR decompositions can fail with errors
        b = aurostd::inverseByQR(a);
      }
      catch (aurostd::xerror& e) {
      }
      if(!aurostd::isequal(a * b, id)) {
        b = aurostd::inverseByAdjoint(a);
      }
      if(!aurostd::isequal(a * b, id)) {
        stringstream message;
        message << "Matrix inversion failed: determinant=" << aurostd::abs(aurostd::det(a)) << " | condition number=" << aurostd::condition_number(a);
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      return b;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>
    bool isNonInvertible(const xmatrix<utype>& a) {  // RETURN ERROR if non invertible  //CO20191201
      utype c=aurostd::det(a);
      if(abs(c)<10e-14) {return TRUE;}
      return FALSE;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  /// @brief calculates inverse of a square matrix using QR decomposition
  ///
  /// @param a square matrix
  ///
  /// @return inverse of the square matrix
  ///
  /// @authors
  /// @mod{SD,20220426,created function}
  template<class utype>                                 // function inverse xmatrix<>
    xmatrix<utype> inverseByQR(const xmatrix<utype>& A) { //SD20220426
      xmatrix<utype> Q, R;
      aurostd::QRDecomposition_HouseHolder(A, Q, R);
      return aurostd::inverseByLUP(R) * trasp(Q);
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  /// @brief calculates LUP decomposition of a square matrix
  ///
  /// @param A square matrix                            
  /// @param LU lower triangular matrix times the upper triangular matrix
  /// @param P permutation matrix
  /// @param tol tolerance when two values are equal
  ///
  /// @authors
  /// @mod{SD,20220426,created function} 
  ///
  /// @note A is a square matrix and LU is the LU decomposition, where LU=(L-I)+U such that A=trasp(P)*LU
  ///
  /// @see
  /// @xlink{https://en.wikipedia.org/wiki/LU_decomposition}
  template<class utype>
    void LUPDecomposition(const xmatrix<utype>& A, xmatrix<double>& LU, xmatrix<double>& P, utype tol) { //SD20220426
      if (!A.issquare) {
        string message = "Matrix needs to be square for LU decomposition";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      LU = aurostd::ones_xm<double>(A.urows, A.ucols, A.lrows, A.lcols);
      for (int i = LU.lrows; i <= LU.urows; i++) {
        for (int j = LU.lcols; j <= LU.ucols; j++) {
          LU(i, j) = (double)A(i, j);
        }
      }
      int imax = 0, _p = 0;
      double maxA = 0.0, absA = 0.0;
      xvector<int> p(LU.urows, LU.lrows);
      xmatrix<double> _LU(LU.urows, LU.ucols, LU.lrows, LU.lcols);
      P = 0.0 * ones_xm<double>(LU.urows, LU.ucols, LU.lrows, LU.lcols);
      for (int i = LU.lrows; i <= LU.urows; i++) {p(i) = i;} // initialize P
      for (int i = LU.lrows; i <= LU.urows; i++) {
        maxA = 0.0;
        imax = i;
        for (int j = i; j <= LU.urows; j++) {
          absA = aurostd::abs(LU(j, i));
          if (absA > maxA) {
            maxA = absA;
            imax = j;
          }
        }
        if (maxA < (double)tol) {
          string message = "Matrix is degenerate";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        }
        if (imax != i) {
          // pivoting P
          _p = p(i);
          p(i) = p(imax);
          p(imax) = _p;
          // pivoting rows of A
          _LU = LU.getxmat(i, i, LU.lcols, LU.ucols);
          LU.setmat(LU.getxmat(imax, imax, LU.lcols, LU.ucols), i, LU.lcols);
          LU.setmat(_LU, imax, LU.lcols);
          //LU.setmat(_LU.getxmat(i, i, LU.lcols, LU.ucols), imax, LU.lcols);
        }
        for (int j = i + 1; j <= LU.urows; j++) {
          LU(j, i) /= LU(i, i);
          for (int k = i + 1; k <= LU.urows; k++) {
            LU(j,k) -= LU(j, i) * LU(i, k);
          }
        }
      }
      for (int i = LU.lrows; i <= LU.urows; i++) {P(i, p(i)) = 1.0;}
    }

  /// @brief calculates LUP decomposition of a square matrix
  ///
  /// @param A square matrix                            
  /// @param L lower triangular matrix
  /// @param U upper triangular matrix
  /// @param P permutation matrix
  /// @param tol tolerance when two values are equal
  ///
  /// @authors
  /// @mod{SD,20220426,created function} 
  ///
  /// @see
  /// @xlink{https://en.wikipedia.org/wiki/LU_decomposition}
  template<class utype>
    void LUPDecomposition(const xmatrix<utype>& A, xmatrix<double>& L, xmatrix<double>& U, xmatrix<double>& P, utype tol) { //SD20220426
      xmatrix<double> LU;
      LUPDecomposition(A, LU, P, tol);
      L = aurostd::eye<double>(LU.urows, LU.ucols, LU.lrows, LU.lcols);
      U = L;
      for (int i = LU.lrows + 1; i <= LU.urows; i++) {
        L.setmat(LU.getxmat(i, i, LU.lcols, i - LU.lcols), i, LU.lcols);
      }
      U += LU - L;
    }

  /// @brief calculates inverse of a square matrix using LUP decomposition
  ///
  /// @param a square matrix
  ///
  /// @return inverse of the square matrix
  ///
  /// @authors
  /// @mod{SD,20220426,created function}
  template<class utype>                                 // function inverse xmatrix<>
    xmatrix<utype> inverseByLUP(const xmatrix<utype>& A) { //SD20220426
      xmatrix<double> LU, P, _IA;
      xmatrix<utype> IA(A.urows, A.ucols, A.lrows, A.lcols);
      LUPDecomposition(A, LU, P);
      _IA = trasp(P);
      for (int i = LU.lrows; i <= LU.urows; i++) {
        for (int j = LU.lcols; j <= LU.ucols; j++) {
          for (int k = 1; k < j; k++) {
            _IA(j, i) -= LU(j, k) * _IA(k, i);
          }
        }
        for (int j = LU.urows; j >= LU.lrows; j--) {
          for (int k = j + 1; k <= LU.urows; k++) {
            _IA(j, i) -= LU(j, k) * _IA(k, i);
          }
          _IA(j, i) /= LU(j, j);
        }
      }
      for (int i = LU.lrows; i <= LU.urows; i++) {
        for (int j = LU.lcols; j <= LU.ucols; j++) {
          IA(i, j) = (utype)_IA(i, j);
        }
      }
      return IA;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  /// @brief calculates condition number of a square matrix
  ///
  /// @param a square matrix
  ///
  /// @return condition number of the matrix
  ///
  /// @authors
  /// @mod{SD,20220426,created function}
  ///
  /// @see
  /// @xlink{https://en.wikipedia.org/wiki/Condition_number}
  template<class utype>
    utype condition_number(const xmatrix<utype>& _a) { //SD20220425
      utype maxv = (utype)0.0, minv = (utype)INFINITY, v;
      xmatrix<utype> a = _a*trasp(_a);
      xvector<utype> rr(a.rows), ri(a.rows);
      aurostd::eigen(a, rr, ri);
      for(int i = 1; i <= a.rows; i++) {
        v = aurostd::modulus(rr(i), ri(i));
        if (v > maxv) {maxv = v;}
        if (v < minv && v > (utype)0.0) {minv = v;}
      }
      return aurostd::sqrt(maxv / minv);
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>  // function roundoff clear small elements
    roundoff(const xmatrix<utype>& a,utype _tol_) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++) {
          if(abs(a[i][j])<(utype) _tol_) c[i][j]=a[i][j]=(utype) 0.0; else c[i][j]=a[i][j];
          //	c[i][j]=nint(a[i][j]/_tol_)*_tol_;
        }
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>  // function roundoff clear small elements
    roundoff(const xmatrix<utype>& a) {
      return roundoff(a,(utype) _AUROSTD_XMATRIX_TOLERANCE_ROUNDOFF_);
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function sum xmatrix<>
    utype modulus(const xmatrix<utype>& a) {  //CO20191110
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function sum: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      return sqrt(modulussquare(a)); //sqrt(sum(trasp(a) * a))
    }
  template<class utype>                                 // function sum xmatrix<>
    utype modulussquare(const xmatrix<utype>& a) {  //CO20191110
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function sum: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      return sum(trasp(a) * a);
    }

  template<class utype>                                 // function sum xmatrix<>
    utype modulus2(const xmatrix<utype>& a) { //CO20191110
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function sum: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      return modulussquare(a);
    }

}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function sum xmatrix<>
    utype sum(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function sum: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      utype c=utype(0.0);
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          c+=a[i][j];
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                           // function sum_colum xmatrix<>
    xvector<utype> sum_column(const xmatrix<utype>& a) {
      xvector<utype> c(a.lcols,a.ucols);
      for(int j=a.lcols;j<=a.ucols;j++) {
        c[j]=(utype) 0.0;
        for(int i=a.lrows;i<=a.urows;i++)
          c[j]+=a[i][j];
      }
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                          // function mean_colum xmatrix<>
    xvector<utype> mean_column(const xmatrix<utype>& a) {
      xvector<utype> c(a.lcols,a.ucols);
      for(int j=a.lcols;j<=a.ucols;j++) {
        c[j]=(utype) 0.0;
        for(int i=a.lrows;i<=a.urows;i++)
          c[j]+=a[i][j];
        c[j]/=(a.urows-a.lrows+1);
      }
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                             // function sum_row xmatrix<>
    xvector<utype> sum_row(const xmatrix<utype>& a) {
      xvector<utype> c(a.lrows,a.urows);
      for(int j=a.lrows;j<=a.urows;j++) {
        c[j]=(utype) 0.0;
        for(int i=a.lcols;i<=a.ucols;i++)
          c[j]+=a[j][i];
      }
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                            // function mean_row xmatrix<>
    xvector<utype> mean_row(const xmatrix<utype>& a) {
      xvector<utype> c(a.lrows,a.urows);
      for(int j=a.lrows;j<=a.urows;j++) {
        c[j]=(utype) 0.0;
        for(int i=a.lcols;i<=a.ucols;i++)
          c[j]+=a[j][i];
        c[j]/=(a.ucols-a.lcols+1);
      }
      return c;
    }
}

// -------------------------------------------------------- functions of xmatrix
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function min xmatrix<>
    utype min(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function min: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      utype c=a[a.lrows][a.lcols];
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          c = c < a[i][j] ? c:a[i][j];
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function min xmatrix<>
    utype min(const xmatrix<utype>& a,int& index_i,int& index_j) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function min: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      utype c=a[a.lrows][a.lcols];
      index_i=a.lrows,index_j=a.lcols;
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          if(a[i][j] < c) {
            c = a[i][j];
            index_i=i;
            index_j=j;
          }
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function max: ");
      printf("index_i=%i, index_j=%i \n",index_i,index_j);
#endif
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function max xmatrix<>
    utype max(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function max: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      utype c=a[a.lrows][a.lcols];
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          c = c > a[i][j] ? c:a[i][j];
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function max xmatrix<>
    utype max(const xmatrix<utype>& a,int& index_i,int& index_j) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function max: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      utype c=a[a.lrows][a.lcols];
      index_i=a.lrows,index_j=a.lcols;
      for(int i=a.lrows;i<=a.urows;i++)
        for(int j=a.lcols;j<=a.ucols;j++)
          if(a[i][j] > c) {
            c = a[i][j];
            index_i=i;
            index_j=j;
          }
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function max: ");
      printf("index_i=%i, index_j=%i \n",index_i,index_j);
#endif
      return c;
    }
}

// ----------------------------------------------------------------------------
//namespace aurostd {
//// namespace aurostd
//template<class utype> double                               // spectral radius
//spectral_radius(const xmatrix<utype>& a)
//{
//#ifdef _XMATH_DEBUG_FUNCTIONS
//printf("M -> function spectral radius: ");
//printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
//printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
//#endif
//if(!a.issquare)
//{cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xmatrix<utype>: failure in spectral radius defined for square xmatrixes" << endl;}
//double out=0.0;
//for(int i=a.lrows;i<=a.urows;i++)
//if(abs(a[i][i])>out)
//out=(double) abs(a[i][i]);
//return out;
//}
//}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> utype  //DX20170115 - double to utype (needed for xcomplex)                                          // trace
    trace(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function trace: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      if(!a.issquare) {
        string message = "Trace is only defined for square matrices";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      utype out=0.0; //DX20170115 - double to utype (needed for xcomplex)
      for(int i=a.lrows;i<=a.urows;i++)
        out+=a[i][i];
      return out;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  /// @brief calculates the Hadamard product of two matrices
  ///
  /// @param A first matrix
  /// @param B second matrix
  ///
  /// @return Hadamard product
  ///
  /// @authors
  /// @mod{SD,20220422,created function} 
  ///
  /// @note Also called element-wise product or Schur product
  ///
  /// @see
  /// @xlink{https://en.wikipedia.org/wiki/Hadamard_product_(matrices)}
  template<class utype>
    xmatrix<utype> HadamardProduct(const xmatrix<utype>& A, const xmatrix<utype>& B) { //SD20220422
      if ((A.rows != B.rows) || (A.cols != B.cols)) {
        stringstream message;
        message << "A and B must have the same dimensions, A.rows=" << A.rows << " B.rows=" << B.rows
                << ", A.cols=" << A.cols << " B.cols=" << B.cols;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      xmatrix<utype> product(A.rows, A.cols);
      for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
          product(i + 1, j + 1) = A(A.lrows + i, A.lcols + j) * B(B.lrows + i, B.lcols + j);
        }
      }
      return product;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {
  template<class utype>
    xmatrix<utype> KroneckerProduct(const xmatrix<utype>& A, const xmatrix<utype>& B) { //ME20180614 - Kronecker product
      int rows, cols, r, c;
      utype aelement, belement;
      rows = A.rows * B.rows;
      cols = A.cols * B.cols;
      xmatrix<utype> product(rows, cols);
      for (int arow = 0; arow < A.rows; arow++) {
        for (int acol = 0; acol < A.cols; acol++) {
          aelement = A(arow + A.lrows, acol + A.lcols);
          for (int brow = 0; brow < B.rows; brow++) {
            for (int bcol = 0; bcol < B.cols; bcol++) {
              belement = B(brow + B.lrows, bcol + B.lcols);
              r = arow * B.rows + 1 + brow;
              c = acol * B.cols + 1 + bcol;
              product(r, c) = aelement * belement;
            }
          }
        }
      }
      return product;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // identity xmatrix
    identity(const xmatrix<utype>& _a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function identity xmatrix: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      if(!_a.issquare) {
        string message = "Identity only defined for square matrices.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      xmatrix<utype> a(_a.urows,_a.ucols,_a.lrows,_a.lcols);
      for(int i=a.lrows;i<=a.urows;i++) a[i][i]=utype(1.0);
      return a;
    }
}

//ME20200123 - if the identity matrix must be square, why would we allow two
// input parameters?
//[OBSOLETE]namespace aurostd {  // namespace aurostd
//[OBSOLETE]  template<class utype>                                 // identity xmatrix
//[OBSOLETE]    xmatrix<utype> identity(const utype& _type,const int& n,const int& m) {
//[OBSOLETE]      xmatrix<utype> a(n,m);
//[OBSOLETE]      if(!a.issquare)
//[OBSOLETE]      {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xmatrix<utype>: failure in identity defined for square xmatrixes [2]" << endl;}
//[OBSOLETE]      for(int i=a.lrows;i<=a.urows;i++)
//[OBSOLETE]        for(int j=a.lcols;j<=a.ucols;j++)
//[OBSOLETE]          if(i==j) a[i][j]=(utype) 1.0; else a[i][j]=(utype) 0.0;
//[OBSOLETE]      return a;
//[OBSOLETE]      if(_type) {;}  // something phony to keep _type busy !
//[OBSOLETE]    }
//[OBSOLETE]}

//ME20200123
namespace aurostd {
  template<class utype>
    xmatrix<utype> identity(const utype& _type, int ubounds, int lbounds) {
      return eye<utype>(ubounds, ubounds, lbounds, lbounds);
      if (_type) {;}  // To suppress compiler warnings
    }
}

//namespace aurostd {  // namespace aurostd
//xmatrix<double>                          // identity_double xmatrix
//identity_double(const int& n,const int& m) {
//xmatrix<double> a(n,m);
//if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xmatrix<utype>: failure in identity_double defined for square xmatrixes [3]" << endl;}
//for(int i=a.lrows;i<=a.urows;i++)
//for(int j=a.lcols;j<=a.ucols;j++)
//if(i==j) a[i][j]=double(1.0); else a[i][j]=double(0.0);
//return a;
//}
//// namespace aurostd
//xmatrix<double>                          // identity_double xmatrix
//identity_double(const int& n) {
//xmatrix<double> a(n,n);
//for(int i=a.lrows;i<=a.urows;i++)
//for(int j=a.lcols;j<=a.ucols;j++)
//if(i==j) a[i][j]=double(1.0); else a[i][j]=double(0.0);
//return a;
//}
//// namespace aurostd
//xmatrix<float>                          // identity_float xmatrix
//identity_float(const int& n,const int& m) {
//xmatrix<float> a(n,m);
//if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xmatrix<utype>: failure in identity_float defined for square xmatrixes [3]" << endl;}
//for(int i=a.lrows;i<=a.urows;i++)
//for(int j=a.lcols;j<=a.ucols;j++)
//if(i==j) a[i][j]=float(1.0); else a[i][j]=float(0.0);
//return a;
//}
//// namespace aurostd
//xmatrix<float>                          // identity_float xmatrix
//identity_float(const int& n) {
//xmatrix<float> a(n,n);
//for(int i=a.lrows;i<=a.urows;i++)
//for(int j=a.lcols;j<=a.ucols;j++)
//if(i==j) a[i][j]=float(1.0); else a[i][j]=float(0.0);
//return a;
//}
//// namespace aurostd
//xmatrix<int>                          // identity_int xmatrix
//identity_int(const int& n,const int& m) {
//xmatrix<int> a(n,m);
//if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xmatrix<utype>: failure in identity_int defined for square xmatrixes [3]" << endl;}
//for(int i=a.lrows;i<=a.urows;i++)
//for(int j=a.lcols;j<=a.ucols;j++)
//if(i==j) a[i][j]=1; else a[i][j]=0;
//return a;
//}
//// namespace aurostd
//xmatrix<int>                          // identity_int xmatrix
//identity_int(const int& n) {
//xmatrix<int> a(n,n);
//for(int i=a.lrows;i<=a.urows;i++)
//for(int j=a.lcols;j<=a.ucols;j++)
//if(i==j) a[i][j]=1; else a[i][j]=0;
//return a;
//}
//// namespace aurostd
//xmatrix<char>                          // identity_char xmatrix
//identity_char(const int& n,const int& m) {
//xmatrix<char> a(n,m);
//if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xmatrix<utype>: failure in identity_char defined for square xmatrixes [3]" << endl;}
//for(int i=a.lrows;i<=a.urows;i++)
//for(int j=a.lcols;j<=a.ucols;j++)
//if(i==j) a[i][j]=1; else a[i][j]=0;
//return a;
//}
//// namespace aurostd
//xmatrix<char>                          // identity_char xmatrix
//identity_char(const int& n) {
//xmatrix<char> a(n,n);
//for(int i=a.lrows;i<=a.urows;i++)
//for(int j=a.lcols;j<=a.ucols;j++)
//if(i==j) a[i][j]=1; else a[i][j]=0;
//return a;
//}
//// namespace aurostd
//xmatrix<bool>                          // identity_bool xmatrix
//identity_bool(const int& n,const int& m) {
//  xmatrix<bool> a(n,m);
//  if(!a.issquare) {cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xmatrix<utype>: failure in identity_bool defined for square xmatrixes [3]" << endl;}
//  for(int i=a.lrows;i<=a.urows;i++)
//    for(int j=a.lcols;j<=a.ucols;j++)
//      if(i==j) a[i][j]=1; else a[i][j]=0;
//  return a;
//}
//// namespace aurostd
//xmatrix<bool>                          // identity_bool xmatrix
//identity_bool(const int& n) {
//  xmatrix<bool> a(n,n);
//  for(int i=a.lrows;i<=a.urows;i++)
//    for(int j=a.lcols;j<=a.ucols;j++)
//      if(i==j) a[i][j]=1; else a[i][j]=0;
//  return a;
//}
//}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  //ME20180904
  template<class utype> xmatrix<utype>                             // conj xmatrix
    conj(const xmatrix<utype>& a){
      if (a.iscomplex) {
        xmatrix<utype> out(a.ucols,a.urows,a.lcols,a.lrows);
        for (int i = a.lrows; i <= a.urows; i++) {
          for (int j = a.lcols; j <= a.ucols; j++) {
            out[i][j] = conj(a[i][j]);
          }
        }
        return out;
      } else {return a;}
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> void
    traspSquareInPlace(xmatrix<utype>& a,bool conjugate){ //CO20191201 - only designed for square matrices
      if(!a.issquare){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::traspSquareInPlace():","this function is designed only for square matrices",_VALUE_ILLEGAL_);}
      int i=0,j=0;
      utype aij=(utype)0;
      if(a.iscomplex&&conjugate){ //CO20191201
        for(i=a.lrows;i<=a.urows;i++){
          for(j=i+1;j<=a.urows;j++){
            aij=a[i][j];
            a[i][j]=conj(a[j][i]);
            a[j][i]=conj(aij);
          }
        }
        return;
      }else{
        for(i=a.lrows;i<=a.urows;i++){
          for(j=i+1;j<=a.urows;j++){
            aij=a[i][j];
            a[i][j]=a[j][i];
            a[j][i]=aij;
          }
        }
        return;
      }
    }
  template<class utype> void
    traspInPlace(const xmatrix<utype>& a,xmatrix<utype>& b,bool conjugate){ //CO20191201
      //always works, just not most efficient for square matrices
      if(b.lrows!=a.lcols || b.lcols!=a.lrows || b.urows!=a.ucols || b.ucols!=a.urows){
        xmatrix<utype> c(a.ucols,a.urows,a.lcols,a.lrows);
        b=c;
      }
      int i=0,j=0;
      if(a.iscomplex&&conjugate){ //CO20191201
        for(i=a.lrows;i<=a.urows;i++){
          for(j=a.lcols;j<=a.ucols;j++){b[j][i]=conj(a[i][j]);}  //cannot do j=i+1 since b is unpopulated
        }
        return;
      }else{
        for(i=a.lrows;i<=a.urows;i++){
          for(j=a.lcols;j<=a.ucols;j++){b[j][i]=a[i][j];}  //cannot do j=i+1 since b is unpopulated
        }
        return;
      }
    }
  template<class utype> void                             // trasp xmatrix
    traspInPlace(xmatrix<utype>& a,bool conjugate){  //ME20190813 - conjugate complex now an option //CO20191201 - in place
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function traspose xmatrix: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      if(a.issquare){traspSquareInPlace(a,conjugate);return;} //this works because we modify a[i][j] and a[j][i] at once
      else{  //cannot trasp in place
        xmatrix<utype> b(a.ucols,a.urows,a.lcols,a.lrows);
        traspInPlace(a,b,conjugate);
        a=b;
        return;
      }
      //[CO20191201 - OBSOLETE]xmatrix<utype> out(a.ucols,a.urows,a.lcols,a.lrows);
      //[CO20191201 - OBSOLETE]for(int i=a.lrows;i<=a.urows;i++)
      //[CO20191201 - OBSOLETE]  for(int j=a.lcols;j<=a.ucols;j++)
      //[CO20191201 - OBSOLETE]    if (a.iscomplex && conjugate) {
      //[CO20191201 - OBSOLETE]      out[j][i] = conj(a[i][j]);
      //[CO20191201 - OBSOLETE]    } else {
      //[CO20191201 - OBSOLETE]      out[j][i]=a[i][j];
      //[CO20191201 - OBSOLETE]    }
      //[CO20191201 - OBSOLETE]return out;
    }
  template<class utype> xmatrix<utype>              // trasp xmatrix
    trasp(const xmatrix<utype>& a,bool conjugate){  //CO20191201
      xmatrix<utype> b(a.ucols,a.urows,a.lcols,a.lrows);
      traspInPlace(a,b,conjugate);
      return b;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                             // trasp xvector
    trasp(const xvector<utype>& a, bool conjugate){  //ME20190813 - conjugate complex now an option
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function traspose xvector: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
#endif
      xmatrix<utype> b(a.urows,1,a.lrows,1);
      int i=0;
      if(a.iscomplex&&conjugate){ //CO20191201
        for(i=a.lrows;i<=a.urows;i++){b[i][1]=conj(a[i]);}
      }else{
        for(i=a.lrows;i<=a.urows;i++){b[i][1]=a[i];}
      }
      return b;
    }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>               // function shift_up xmatrix<>
    shift_up(const xmatrix<utype>& a) {
      utype aus;
      for(int j=a.lcols;j<=a.ucols;j++) {
        aus=a[a.lrows][j];
        for(int i=a.lrows;i<=a.urows-1;i++)
          a[i][j]=a[i+1][j];
        a[a.urows][j]=aus;
      }
      return a;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>             // function shift_down xmatrix<>
    shift_down(const xmatrix<utype>& a) {
      utype aus;
      for(int j=a.lcols;j<=a.ucols;j++) {
        aus=a[a.urows][j];
        for(int i=a.urows;i>=a.lrows+1;i--)
          a[i][j]=a[i-1][j];
        a[a.lrows][j]=aus;
      }
      return a;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>             // function shift_left xmatrix<>
    shift_left(const xmatrix<utype>& a) {
      utype aus;
      for(int i=a.lrows;i<=a.urows;i++) {
        aus=a[i][a.lcols];
        for(int j=a.lcols;j<=a.ucols-1;j++)
          a[i][j]=a[i][j+1];
        a[i][a.ucols]=aus;
      }
      return a;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>             // function shift_right xmatrix<>
    shift_right(const xmatrix<utype>& a) {
      utype aus;
      for(int i=a.lrows;i<=a.urows;i++) {
        aus=a[i][a.ucols];
        for(int j=a.ucols;j>=a.lcols+1;j++)
          a[i][j]=a[i][j+1];
        a[i][a.lcols]=aus;
      }
      return a;
    }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
// SWAP THINGS
namespace aurostd {  // namespace aurostd
  template<class utype> void                                // swap_columns
    swap_cols(xmatrix<utype>& a,const int& i,const int& j) {  // swap_columns
      if(i<a.lcols || i>a.ucols) return; // nothing to do, out of boundaries
      if(j<a.lcols || j>a.ucols) return; // nothing to do, out of boundaries
      if(i==j) return; // nothing to do, no swap
      utype temp;
      for(int k=a.lrows;k<=a.urows;k++) {
        temp=a[k][i];a[k][i]=a[k][j];a[k][j]=temp;
      }
    }

  template<class utype> void                                  // swap_columns
    swap_columns(xmatrix<utype>& a,const int& i,const int& j) {  // swap_columns
      swap_cols(a,i,j);
      //     if(i<a.lcols || i>a.ucols) return; // nothing to do, out of boundaries
      //     if(j<a.lcols || j>a.ucols) return; // nothing to do, out of boundaries
      //     if(i==j) return; // nothing to do, no swap
      //     utype temp;
      //     for(int k=a.lrows;k<=a.urows;k++) {
      //       temp=a[k][i];a[k][i]=a[k][j];a[k][j]=temp;
      //     }
    }

  template<class utype> void                                // swap_rows
    swap_rows(xmatrix<utype>& a,const int& i,const int& j) {  // swap_rows
      if(i<a.lrows || i>a.urows) return; // nothing to do, out of boundaries
      if(j<a.lrows || j>a.urows) return; // nothing to do, out of boundaries
      if(i==j) return; // nothing to do, no swap
      utype temp;
      for(int k=a.lcols;k<=a.ucols;k++) {
        temp=a[i][k];a[i][k]=a[j][k];a[j][k]=temp;
      }
    }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // sign xmatrix
    sign(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=aurostd::sign(a[i][j]);
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // nint xmatrix
    nint(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=aurostd::nint(a[i][j]);
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // floor xmatrix
    floor(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=std::floor(a[i][j]);
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // trunc xmatrix
    trunc(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          //	c[i][j]=std::trunc(a[i][j]);
          c[i][j]=trunc(a[i][j]);
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // round xmatrix
    round(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          //	c[i][j]=std::round(a[i][j]);
          c[i][j]=round(a[i][j]);
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                          // ceil xmatrix
    ceil(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=std::ceil(a[i][j]);
      return c;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                           // mabs xmatrix
    mabs(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=aurostd::abs(a[i][j]);
      return c;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                           // abs xmatrix
    abs(const xmatrix<utype>& a) {
      xmatrix<utype> c(a.urows,a.ucols,a.lrows,a.lcols);
      for(int i=c.lrows;i<=c.urows;i++)
        for(int j=c.lcols;j<=c.ucols;j++)
          c[i][j]=aurostd::abs(a[i][j]);
      return c;
    }
}

// ****************************************************************************
// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                               // exp xmatrix
    exp_old(const xmatrix<utype>& a) {
      if(!a.issquare) {
        string message = "exp only defined for square matrices.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
      // UNUSED   bool convergence=FALSE;
      for(int n=0;n<=1000;n++) {
        if(n==0) identity(an);
        else an=(an*a)/utype(n);
        out=out+an;
        //    if(abs(trace(an)/trace(out))<_exponential_convergence)
      }
      return out;
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                               // exp xmatrix
    exp(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function exponential xmatrix: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      if(!a.issquare) {
        string message = "exp only defined for square matrices.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
      bool convergence=FALSE;
      for(int n=0;!convergence;n++) {
        if(n==0) identity(an);
        else an=(an*a)/utype(n);
        out=out+an;
        // cerr << n << endl;
        // if(abs(trace(an)/trace(out))<_exponential_convergence)
        if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
        if(n>100) convergence=TRUE;
      }
      return out;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                               // sin xmatrix
    sin(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function sin xmatrix: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      if(!a.issquare) {
        string message = "sin only defined for square matrices.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols), an(a.urows,a.ucols,a.lrows,a.lcols);
      bool convergence=FALSE;
      for(int n=0;!convergence;n++) {
        if(n==0) an=a;
        else an=(a*a*an)/utype(-1.0*(2.0*n+1.0)*(2.0*n));
        out=out+an;
        //    if(abs(trace(an)/trace(out))<_exponential_convergence)
        if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
        if(n>100) convergence=TRUE;
      }
      return out;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                               // cos xmatrix
    cos(const xmatrix<utype>& a) {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function cos xmatrix: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      if(!a.issquare) {
        string message = "cos only defined for square matrices.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
      bool convergence=FALSE;
      for(int n=0;!convergence;n++) {
        if(n==0) identity(an);
        else an=(a*a*an)/utype(-1.0*(2.0*n)*(2.0*n-1.0));
        out=out+an;
        if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
        if(n>100) convergence=TRUE;
      }
      return out;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                              // sinh xmatrix
    sinh(const xmatrix<utype>& a)
    {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function sinh xmatrix: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      if(!a.issquare) {
        string message = "sinh only defined for square matrices.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
      bool convergence=FALSE;
      for(int n=0;!convergence;n++) {
        if(n==0) an=a;
        else an=(a*a*an)/utype((2.0*n+1.0)*(2.0*n));
        out=out+an;
        // if(abs(trace(an)/trace(out))<_exponential_convergence)
        if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
        if(n>100) convergence=TRUE;
      }
      return out;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xmatrix<utype>                              // cosh xmatrix
    cosh(const xmatrix<utype>& a)
    {
#ifdef _XMATH_DEBUG_FUNCTIONS
      printf("M -> function cosh xmatrix: ");
      printf("a.lrows=%i, a.urows=%i, ",a.lrows,a.urows);
      printf("a.lcols=%i, a.ucols=%i\n",a.lcols,a.ucols);
#endif
      if(!a.issquare) {
        string message = "cosh only defined for square matrices.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      xmatrix<utype> out(a.urows,a.ucols,a.lrows,a.lcols),an(a.urows,a.ucols,a.lrows,a.lcols);
      bool convergence=FALSE;
      for(int n=0;!convergence;n++) {
        if(n==0) identity(an);
        else an=(a*a*an)/utype((2.0*n)*(2.0*n-1.0));
        out=out+an;
        // if(abs(trace(an)/trace(out))<_exponential_convergence)
        if(n>30) if(abs(trace(an))<_exponential_convergence) convergence=TRUE;
        if(n>100) convergence=TRUE;
      }
      return out;
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  /// @brief solves A*x=b, where A is a m-by-n matrix and b is a m-by-1 vector and
  /// the output x is a n-by-1 vector
  ///
  /// @param A m-by-n matrix
  /// @param b m-by-1 vector
  ///
  /// @return least squares approximation as a n-by-1 vector
  ///
  /// @authors
  /// @mod{SD,20220426,created function}
  ///
  /// @see
  /// @xlink{https://en.wikipedia.org/wiki/Linear_least_squares}
  template<class utype>                                   
    xvector<utype> LinearLeastSquares(const xmatrix<utype>& A, const xvector<utype>& b) {
    return aurostd::inverse(trasp(A) * A) * trasp(A) * b;
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // GaussJordan xmatrix
    void GaussJordan(xmatrix<utype>& A, xmatrix<utype>& B) {
      /// This function uses Gaussian Jordan elimination to solve A*x=b.  It returns the solution x and the inverse of A.
      string message = "";
      if(A.lrows!=1) {
        message = "[1] A.lrows!=1 <<  A.lrows=" + aurostd::utype2string<int>(A.lrows);
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_ERROR_);
      }
      if(A.lcols!=1) {
        message = "[2] A.lcols!=1 <<  A.lcols=" + aurostd::utype2string<int>(A.lcols);
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_ERROR_);
      }
      if(B.lrows!=1) {
        message = "[3] B.lrows!=1 <<  B.lrows=" + aurostd::utype2string<int>(B.lrows);
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_ERROR_);
      }
      if(B.lcols!=1) {
        message = "[4] B.lcols!=1 <<  B.lcols=" + aurostd::utype2string<int>(B.lcols);
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_ERROR_);
      }
      if(A.urows!=A.ucols) {
        message = "[5] A.urows!=A.ucols <<  A.urows=" + aurostd::utype2string<int>(A.urows) + " A.ucols=" + aurostd::utype2string<int>(A.ucols);
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      if(A.ucols!=B.urows) {
        message = "[6] A.ucols!=B.urows <<  A.ucols=" + aurostd::utype2string<int>(A.ucols) + " B.urows=" + aurostd::utype2string<int>(B.urows);
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      int n=A.urows;
      int m=B.ucols;

      // cerr << "GaussJordan" << A.urows << " " << A.ucols << endl;

      int i,icol=1,irow=1,j,k,l,ll;
      utype big,dum,pivinv,temp;

      xvector<int> indxc(n),indxr(n),ipiv(n);

      for(j=1;j<=n;j++) ipiv[j]=0;
      for(i=1;i<=n;i++) {
        big=0.0;
        for(j=1;j<=n;j++)
          if(ipiv[j]!=1)
            for(k=1;k<=n;k++) {
              if(ipiv[k] == 0) {
                if(aurostd::abs(A[j][k])>=big) {
                  big=aurostd::abs(A[j][k]);
                  irow=j;
                  icol=k;
                }
              } else if(ipiv[k]>1) {
                message = "[7]: Singular Matrix-1";
                throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
              }
            }
        ++(ipiv[icol]);
        if(irow!=icol) {
          for(l=1;l<=n;l++) SWAP(A[irow][l],A[icol][l]);
          for(l=1;l<=m;l++) SWAP(B[irow][l],B[icol][l]);
        }
        indxr[i]=irow;
        indxc[i]=icol;
        if(A[icol][icol]==(double) 0.0) {
          message = "[8]: Singular Matrix-2";
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
        pivinv=1.0/A[icol][icol];
        A[icol][icol]=1.0;
        for(l=1;l<=n;l++) A[icol][l]*=pivinv;
        for(l=1;l<=m;l++) B[icol][l]*=pivinv;
        for(ll=1;ll<=n;ll++)
          if(ll!=icol) {
            dum=A[ll][icol];
            A[ll][icol]=0.0;
            for(l=1;l<=n;l++) A[ll][l]-=A[icol][l]*dum;
            for(l=1;l<=m;l++) B[ll][l]-=B[icol][l]*dum;
          }
      }
      for(l=n;l>=1;l--) {
        if(indxr[l]!=indxc[l])
          for(k=1;k<=n;k++)
            SWAP(A[k][indxr[l]],A[k][indxc[l]]);
      }
    }
}

// ----------------------------------------------------------------------------
namespace aurostd {   // least square stuff aurostd adaptation of nrecipes    // 1 August 2014
  // namespace aurostd
  template<class utype> void gaussj(xmatrix<utype>& a, int n, xmatrix<utype>& b, int m) {  // with indices
    // linear equation solution by gauss-jordan elimination, a[1,n][1,n] is the input matrix.
    // b[1,n][1,m] is input containing the m right-hand side vectors. On the output a is replaced
    // by its matrix inverse, and b is replaced by the corresponding set of solution vectors.
    int i,icol=0,irow=0,j,k,l,ll; // default definitions to avoid compilation errors
    utype big,dum,pivinv,temp;

    string message = "";
    if(n>a.rows) {
      message = "n>a.rows";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_RANGE_);
    }
    if(n>b.rows) {
      message = "n>b.rows";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_RANGE_);
    }
    if(m>b.cols) {
      message = "m>b.cols";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_RANGE_);
    }

    xvector<int> indxc(1,n);
    xvector<int> indxr(1,n);
    xvector<int> ipiv(1,n);
    for (j=1;j<=(int) n;j++) ipiv[j]=0;
    for (i=1;i<=(int) n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
        if(ipiv[j] != 1)
          for (k=1;k<=n;k++) {
            if(ipiv[k] == 0) {
              if(aurostd::abs(a[j][k]) >= big) {
                big=aurostd::abs(a[j][k]);
                irow=j;
                icol=k;
              }
            } else if(ipiv[k] > 1) {
              message = "Singular Matrix-1";
              throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
            }
          }
      ++(ipiv[icol]);
      if(irow != icol) {
        for (l=1;l<=n;l++) {SWAP(a[irow][l],a[icol][l]);}
        for (l=1;l<=m;l++) {SWAP(b[irow][l],b[icol][l]);}
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if(a[icol][icol] == 0.0) {
        message = "Singular Matrix-2";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=1;l<=n;l++) a[icol][l] *= pivinv;
      for (l=1;l<=m;l++) b[icol][l] *= pivinv;
      for (ll=1;ll<=n;ll++)
        if(ll != icol) {
          dum=a[ll][icol];
          a[ll][icol]=0.0;
          for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
          for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
        }
    }
    for (l=n;l>=1;l--) {
      if(indxr[l] != indxc[l])
        for (k=1;k<=n;k++) {
          SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
    }
  }
}

// ----------------------------------------------------------------------------
namespace aurostd {   // least square stuff aurostd adaptation of nrecipes    // 1 August 2014
  template<class utype>
    void lfit(xvector<utype> x, xvector<utype> y, xvector<utype> sig,
        xvector<utype>& a, xvector<int> ia,
        xmatrix<utype>& covar, utype& chisq,
        void (*funcs)(utype, xvector<utype>&)) {
      // Given a set of data points x[1,ndat],y[1,ndat] with individual standar deviation sig[1,ndat], use chisq minimization to fit for some or all the coefficients
      // a[1,ma] of a function that depends linearly on a, y=sum_i a_i*afunc_i(x). The input array ia[1,ma] indicates by nonzero entries those componends of a
      // that should be fitted for, and by zero entries those components that should be held fiuxed at their input values.
      // The prgram returns value for a[1,ma], chisq,  and the covariance atrix covar[1,ma][1,ma]. (Parameters held fixed will return zero covariances.).
      // The user supplies a routine funcs(x,xvector<afunc>) that returns the ma basis funcions evaluated at x=X in the array afunc[1,ma]

      string message = "";
      int ndat=x.rows;
      if(y.rows!=x.rows) {
        message = "y.rows!=x.rows";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      if(sig.rows!=x.rows) {
        message = "sig.rows!=x.rows";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }

      int ma=a.rows;
      if(ia.rows!=a.rows) {
        message = "ia.rows!=a.rows";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }

      int i,j,k,l,m,mfit=0;
      utype ym,wt,sum,sig2i;

      aurostd::xmatrix<utype> beta(1,ma,1,1);
      aurostd::xvector<utype> afunc(1,ma);
      for (j=1;j<=(int) ma;j++)
        if(ia[j]) mfit++;
      if(mfit == 0) {
        message = "no parameters to be fitted";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      for (j=1;j<=mfit;j++) {
        for (k=1;k<=mfit;k++) covar[j][k]=0.0;
        beta[j][1]=0.0;
      }
      for (i=1;i<=(int) ndat;i++) {
        (*funcs)(x[i],afunc);
        ym=y[i];
        if(mfit < (int) ma) {
          for (j=1;j<=(int) ma;j++)
            if(!ia[j]) ym -= a[j]*afunc[j];
        }
        sig2i=1.0/(sig[i]*sig[i]);
        for (j=0,l=1;l<=(int) ma;l++) {
          if(ia[l]) {
            wt=afunc[l]*sig2i;
            for (j++,k=0,m=1;m<=l;m++)
              if(ia[m]) covar[j][++k] += wt*afunc[m];
            beta[j][1] += ym*wt;
          }
        }
      }
      for (j=2;j<=mfit;j++)
        for (k=1;k<j;k++)
          covar[k][j]=covar[j][k];
      gaussj(covar,mfit,beta,1); // operate up to mfit
      for (j=0,l=1;l<=(int) ma;l++)
        if(ia[l]) a[l]=beta[++j][1];
      chisq=0.0;
      for (i=1;i<=(int) ndat;i++) {
        (*funcs)(x[i],afunc);
        for (sum=0.0,j=1;j<=(int) ma;j++) sum += a[j]*afunc[j];
        chisq += (((y[i]-sum)/sig[i])*((y[i]-sum)/sig[i]));
      }
      covsrt(covar,ia,mfit);
    }

}
// ----------------------------------------------------------------------------
namespace aurostd {   // least square stuff aurostd adaptation of nrecipes    // 1 August 2014
  template<class utype> void covsrt(xmatrix<utype>&covar, xvector<int> ia, int mfit) {
    // Expand in storage the covariance matrix covar[1,ma][1,ma], so as to take into account parameters
    // that are being fixed. (for the latter, return zero covariance.)
    int i,j,k;
    utype temp;

    int ma=covar.rows;
    if(covar.cols!=covar.rows) {
      string message = "covar.cols!=covar.rows";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    for (i=mfit+1;i<=ma;i++) {
      for (j=1;j<=i;j++) {
        covar[i][j]=covar[j][i]=0.0;
      }
    }
    k=mfit;
    for (j=ma;j>=1;j--) {
      if(ia[j]) {
        for (i=1;i<=ma;i++) {SWAP(covar[i][k],covar[i][j]);}
        for (i=1;i<=ma;i++) {SWAP(covar[k][i],covar[j][i]);}
        k--;
      }
    }
  }
}

namespace aurostd {  // namespace aurostd
  void GCD(const xmatrix<int>& ma,const xmatrix<int>& mb,xmatrix<int>& mgcd){ //CO20191201
    if(ma.rows==0 || ma.cols==0){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::GCD():","ma.rows==0 || ma.cols==0",_INPUT_NUMBER_);}
    xmatrix<int> mx(ma.urows,ma.ucols,ma.lrows,ma.lcols),my(ma.urows,ma.ucols,ma.lrows,ma.lcols);
    return GCD(ma,mb,mgcd,mx,my);
  }
  void GCD(const xmatrix<int>& ma,const xmatrix<int>& mb,xmatrix<int>& mgcd,xmatrix<int>& mx,xmatrix<int>& my){ //CO20191219
    if(ma.rows==0 || ma.cols==0){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::GCD():","ma.rows==0 || ma.cols==0",_INPUT_NUMBER_);}
    //ma vs. mb
    if(ma.lrows!=mb.lrows){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::GCD():","ma.lrows!=mb.lrows",_INDEX_MISMATCH_);}
    if(ma.urows!=mb.urows){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::GCD():","ma.urows!=mb.urows",_INDEX_MISMATCH_);}
    if(ma.lcols!=mb.lcols){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::GCD():","ma.lcols!=mb.lcols",_INDEX_MISMATCH_);}
    if(ma.ucols!=mb.ucols){throw aurostd::xerror(__AFLOW_FILE__,"aurostd::GCD():","ma.ucols!=mb.ucols",_INDEX_MISMATCH_);}
    //ma vs. mgcd
    if(ma.lrows!=mgcd.lrows || ma.urows!=mgcd.urows || ma.lcols!=mgcd.lcols || ma.ucols!=mgcd.ucols){xmatrix<int> mgcd_tmp(ma);mgcd=mgcd_tmp;}
    //ma vs. mx
    if(ma.lrows!=mx.lrows || ma.urows!=mx.urows || ma.lcols!=mx.lcols || ma.ucols!=mx.ucols){xmatrix<int> mx_tmp(ma);mx=mx_tmp;}
    //ma vs. my
    if(ma.lrows!=my.lrows || ma.urows!=my.urows || ma.lcols!=my.lcols || ma.ucols!=my.ucols){xmatrix<int> my_tmp(ma);my=my_tmp;}
    for(int i=ma.lrows;i<=ma.urows;i++){
      for(int j=ma.lcols;j<=ma.ucols;j++){
        GCD(ma[i][j],mb[i][j],mgcd[i][j],mx[i][j],my[i][j]);
      }
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// DX20201125
namespace aurostd {
  template<class utype> void polarDecomposition(const xmatrix<utype>& transformation_matrix,
      xmatrix<utype>& rotation,
      xmatrix<utype>& deformation,
      bool check_orthogonal_rotation) {

    // Decompose a transformation into its rotation and deformation
    // matrices: T=R*U (where T=original matrix, R=rotation, U=deformation).
    // Procedure:
    // 1) T^2 = trasp(T)*T=trasp(R*U)*(R*U)=trasp(U)*trasp(R)*R*U=trasp(U)*U
    //    (since trasp(R)*R=I, i.e. orthogonal matrix)
    // 2) U = sqrt(trasp(U)*U), using diagonalization technique:
    //    http://en.wikipedia.org/wiki/Square_root_of_a_matrix#By_diagonalization
    // 3) R = T*inverse(U)
    // Following ref: http://www.continuummechanics.org/polardecomposition.html
    // Generalized for an nxn matrix.

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;

    // ---------------------------------------------------------------------------
    // analysis only works for a square matrix
    if(!transformation_matrix.issquare){
      message << "The transformation matrix must be a square matrix.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
    }

    // ---------------------------------------------------------------------------
    // save matrix dimension
    uint dimension = transformation_matrix.rows;

    // ---------------------------------------------------------------------------
    // square the matrix
    xmatrix<utype> T_squared = trasp(transformation_matrix)*transformation_matrix;

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " T^2: " << T_squared << endl; }

    // ---------------------------------------------------------------------------
    // if T^2 is the identity, then this is a unitary transformation, no deformation
    // (not sure how sensitive the deformation is, we may not be able to do this)
    if(aurostd::isidentity(T_squared)){
      rotation = transformation_matrix;
      deformation = aurostd::eye<utype>(dimension,dimension); //DX+ME20210111
      return;
    }

    // ---------------------------------------------------------------------------
    // find square root of matrix via diagonalization method
    xvector<utype> diag(dimension); //diag: vector of diagonal components
    xmatrix<utype> eigen_vec(dimension,dimension); //eigen_vec: matrix with eigen vectors as columns

    // ---------------------------------------------------------------------------
    // Jacobi
    jacobi(T_squared, diag, eigen_vec);

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " diag: " << diag << endl;
      cerr << __AFLOW_FUNC__ << " eigen_vec: " << eigen_vec << endl;
    }

    // ---------------------------------------------------------------------------
    // build diagonal matrix
    xmatrix<utype> diag_matrix = aurostd::eye<utype>(dimension,dimension);
    for(uint i=1;i<=dimension;i++){
      diag_matrix[i][i] = aurostd::sqrt(diag(i));
    }

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " diag_matrix: " << diag_matrix << endl; }

    // ---------------------------------------------------------------------------
    // find deformation (U) via U=v*D*inverse(v);
    deformation = eigen_vec*diag_matrix*aurostd::inverse(eigen_vec);

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " deformation matrix (U): " << deformation << endl; }

    // ---------------------------------------------------------------------------
    // find rotation (R) via R=T*inverse(U)
    rotation = transformation_matrix*aurostd::inverse(deformation);

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " rotation matrix (R): " << rotation << endl; }

    // ---------------------------------------------------------------------------
    // verify conditions of an orthogonal matrix for rotation
    if(check_orthogonal_rotation){
      xmatrix<utype> identity_matrix = rotation*trasp(rotation);
      if(LDEBUG){
        // R^T==R^-1
        cerr << __AFLOW_FUNC__ << " transpose(R):" << endl << aurostd::trasp(rotation) << endl;
        cerr << __AFLOW_FUNC__ << " inverse(R):" << endl << aurostd::inverse(rotation) << endl;
        // R*R^T=I
        cerr << __AFLOW_FUNC__ << " identity? (R*R^T=I):" << endl << identity_matrix << endl;
      }
      if(!aurostd::isidentity(identity_matrix)){
        message << "Extracted rotation should be an orthogonal matrix (R*R^T==I):" << endl << identity_matrix;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
  }
}


// ----------------------------------------------------------------------------
namespace aurostd {
  /// @brief equilibrates a matrix such that it has max-norm 1 in every row and column
  ///
  /// @param A_orig original matrix
  /// @param A equilibrated matrix
  /// @param R row scaling matrix
  /// @param C column scaling matrix
  /// @param tol tolerance when the equilibration is converged
  ///
  /// @authors
  /// @mod{SD,20220425,created function}
  ///
  /// @note The equilibrated matrix A can be written in terms of the original matrix as: A = R * A_orig * C
  ///
  /// @see
  /// @doi{10.1137/S0895479891222088}
  /// @xlink{https://cs.stanford.edu/people/paulliu/files/cs517-project.pdf}
  template<class utype> void equilibrateMatrix(const xmatrix<utype>& A_orig, xmatrix<utype>& A, xmatrix<utype>& R, xmatrix<utype>& C, uint niter, utype tol) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    A = A_orig;
    A.shift(1, 1);
    int m = A.rows, n = A.cols, mn = m + n;
    //           ( 0 | A )
    // A_block = ( ----- ) // symmetric matrix
    //           ( A'| 0 )
    xmatrix<utype> A_block = (utype)0.0 * aurostd::ones_xm<utype>(mn, mn);
    A_block.setmat(A, 1, m + 1);
    A_block.setmat(trasp(A), m + 1, 1);
    A = A_block;
    xmatrix<utype> D = aurostd::eye<utype>(mn), D_tmp;
    utype d;
    bool flag_conv = false;
    uint iter = 0;
    if (LDEBUG) {cerr << __AFLOW_FUNC__ << " A_orig=" << A_orig << endl;}
    while (iter < niter && flag_conv != true) {
      flag_conv = true;
      for (int i = 1; i <= mn; i++) {
        D_tmp = aurostd::eye<utype>(mn);
        d = (utype)1.0 / aurostd::sqrt(aurostd::linf_norm(A.getxmat(1, mn, i, i)));
        D_tmp(i, i) = d;
        D = D * D_tmp;
        A = D_tmp * A * D_tmp;
        if (flag_conv && (aurostd::abs(1.0 - d) > tol)) {flag_conv = false;}
      }
      iter++;
    }
    R = D.getxmat(1, m, 1, m);
    C = D.getxmat(m + 1, mn, m + 1, mn);
    A = R * A_orig * C;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " iter=" << iter << endl;
      cerr << __AFLOW_FUNC__ << " A=" << endl << A << endl;
    }
  }
}

// ----------------------------------------------------------------------------

namespace aurostd {  // namespace aurostd
  template<class utype>                                      //  std::cout operator <<
    std::ostream& operator<< (std::ostream& buf,const xmatrix<utype>& x) {
      char buf2[80];
      string iobuf1,iobuf2,iobuf3,iobuf4,iobuf5;         // buffers
      utype xij=0;int i,j;                                             // buffer x[i]
      bool done=FALSE;
      if(_isfloat(xij)) {                                    // floating point mode
        if(!_iscomplex(xij)) {                                     // real numbers
          if(_size(xij)==sizeof(long double)) {                     // long double
            iobuf1="%13.7lle";                                     // long double
            iobuf2="long double  ";                                  // long double
            iobuf3="    0.0      ";                                // long double
            iobuf4="[%2i]   ";                                       // long double
            iobuf5="     ";                                        // long double
            done=TRUE;
          }	
          if(_size(xij)==sizeof(double)) {                               // double
            iobuf1="%11.4le";                                           // double
            iobuf2="double     ";                                         // double
            iobuf3="   0.0     ";                                       // double
            iobuf4="[%2i]  ";                                             // double
            iobuf5="    ";                                              // double
            done=TRUE;
          }	
          if(_size(xij)==sizeof(float)) {                                 // float
            iobuf1=" % 2.4lf";                                             // float
            iobuf2="float    ";                                            // float
            iobuf3="  0.0000";                                             // float
            iobuf4="[%2i] ";                                               // float
            iobuf5="    ";
            done=TRUE;
          }	
        } else {                                                // xcomplex numbers
          if(_size(xij)==sizeof(xcomplex<long double>)) {  // xcomplex<long double>
            // 	  iobuf1=" (% 13.7lle,% 13.7lle)";                // xcomplex<long double>
            // 	  iobuf2="xcomplex<long double>";                 // xcomplex<long double>
            // 	  iobuf3=" (     0.0      ,     0.0      )";      // xcomplex<long double>
            // 	  iobuf4="[%2i]   ";                              // xcomplex<long double>
            // 	  iobuf5="       ";                               // xcomplex<long double>
            done=TRUE;
          }	
          if(_size(xij)==sizeof(xcomplex<double>)) {            // xcomplex<double>
            // 	  iobuf1=" (% 11.5le,% 11.5le)";                       // xcomplex<double>
            // 	  iobuf2="xcomplex<double>";                           // xcomplex<double>
            // 	  iobuf3=" (   0.0      ,   0.0      )";               // xcomplex<double>
            // 	  iobuf4="[%2i]  ";                                    // xcomplex<double>
            // 	  iobuf5="      ";                                     // xcomplex<double>
            done=TRUE;
          }	
          if(_size(xij)==sizeof(xcomplex<float>)) {              // xcomplex<float>
            // 	  iobuf1=" (% 10.4le,% 10.4le)";                        // xcomplex<float>
            // 	  iobuf2="xcomplex<float>";                             // xcomplex<float>
            // 	  iobuf3=" (   0.0     ,   0.0     )";                  // xcomplex<float>
            // 	  iobuf4="[%2i] ";                                      // xcomplex<float>
            // 	  iobuf5="      ";                                      // xcomplex<float>
            done=TRUE;
          }	
        }
      } else {                                                      // integer mode
        if(_size(xij)==sizeof(long int))  {                            // long int
          iobuf1="%11i";                                                // long int
          iobuf2="long int     ";                                       // long int
          iobuf3="        0 ";                                          // long int
          iobuf4="[%2i] ";                                              // long int
          iobuf5="     ";                                               // long int
          done=TRUE;
        }	
        if(_size(xij)==sizeof(int))  {                                      // int
          iobuf1="%11i";                                                     // int
          iobuf2="int          ";                                            // int
          iobuf3="        0 ";                                               // int
          iobuf4="[%2i] ";                                                   // int
          iobuf5="     ";                                                    // int
          done=TRUE;
        }	
        if(_size(xij)==sizeof(char))  {                                    // char
          iobuf1="%3d";                                                     // char
          iobuf2="char ";                                                   // char
          iobuf3="  0 ";                                                    // char
          iobuf4="[%2i] ";                                                  // char
          iobuf5="";                                                        // char
          done=TRUE;
        }	
      }

      //cerr << iobuf2 << endl;
      if(done==FALSE) {
        string message = "no data type available for user type";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

#ifdef _XMATH_DEBUG_OUTPUT
      buf << iobuf2;
      for(j=x.lcols;j<=x.ucols;j++) {
        snprintf(buf2,80,iobuf4.c_str(),j);                                            // above
        buf << buf2 << iobuf5;
      }
      buf << endl ;
#endif
      for(i=x.lrows;i<=x.urows;i++) {
#ifdef _XMATH_DEBUG_OUTPUT
        snprintf(buf2,80,iobuf4.c_str(),i);                                             // near
        buf << buf2 ;
#endif
#ifdef _XMATH_LATGEN_AL_GULP
        buf << "Al core" ;
#endif
        for(j=x.lcols;j<=x.ucols;j++) {
          xij=x[i][j];
          if(!_iscomplex(xij)) {
            if(_isfloat(xij)) {
              if(abs(xij)> (double) _xmatrix_epsilon) {
                snprintf(buf2,80,iobuf1.c_str(),aurostd::_real(xij));
              } else {
                //	      sprintf(buf2,iobuf3);
                snprintf(buf2,80,iobuf1.c_str(),aurostd::_real(xij));
              }
            } else {
              if(aurostd::_real(xij)!=0) {
                snprintf(buf2,80,iobuf1.c_str(),aurostd::_real(xij));
              } else {
                //  sprintf(buf2,iobuf3.c_str());
                snprintf(buf2,80,iobuf1.c_str(),aurostd::_real(xij));
              }
            }
            buf << string(buf2) << " ";
          } else {
            //	  if(abs(xij)>  (float) _xmatrix_epsilon)
            //    sprintf(buf2,iobuf1.c_str(),real(xij),imag(xij));
            //  else
            //    sprintf(buf2,iobuf3.c_str());
            buf << xij << " ";  // problem of printing in xcomplex
          }
          //	if(j<x.ucols) buf;
        }
        if(i<x.urows) buf << endl;
      }
      // cerr << "[3]" << endl;
      return buf;
    }
}


// ----------------------------------------------------------------------------
// --------------------------------------------- MODULE vector_matrix_utilities


//*****************************************************************************
// ORTHOGONALITY STUFF
#define _DEFAULT_EPS_BASIS_REDUCE_ 0.001

namespace aurostd {
  // namespace aurostd

  // ***************************************************************************************************
  //  This function calculates the "orthogonality defect" of the given basis of a 3D lattice.
  template<class utype> utype
    orthogonality_defect(const xmatrix<utype>& basis) {
      utype od=1.0;
      xvector<utype> bj(3);
      for(int j=1;j<=3;j++) {
        bj(1)=basis(1,j);
        bj(2)=basis(2,j);
        bj(3)=basis(3,j);
        od=od*aurostd::modulus(bj);
      }
      od=od/aurostd::abs(aurostd::det(basis));
      return od;
    }

  // ***************************************************************************************************
  // This routine takes two vectors (in three-space) and reduces them to form a shortest set (Minkowski
  // reduced). The idea is to subtract B from C so that the new C is as close to the origin as any
  // lattice point along the line that passes through C in the direction of B. The process is repeated
  // then for C subtracted from B, and so on, until the new vector isn't shorter than the other. It's
  // pretty obvious if you do an example by hand. Also see 3.1 of Lecture notes in computer science,
  // ISSN 0302-974, ANTS - VI : algorithmic number theory, 2004, vol. 3076, pp. 338-357 ISBN
  // 3-540-22156-5
  template<class utype> bool
    gaussian_reduce_two_vectors(xvector<utype>& B, xvector<utype>& C,utype eps) {
      xvector<utype> temp(3);
      for(int it=0;;it++) { // dont touch this
        int SwapCnt=0;//GH Counter for the number of times B and C are swapped
        if(it > 100) { cerr << "gaussian_reduce_two_vectors failed to converge" << endl;return FALSE;}
        if(aurostd::modulus(B) > aurostd::modulus(C)) {
          temp=B;  // Keep C as the longest vector
          B=C;     // Keep C as the longest vector
          C=temp;  // Keep C as the longest vector
          SwapCnt++; //GH Keep track of the number of swaps
        }
        //    cerr << aurostd::modulus(C) << " " << scalar_product(B,C)/aurostd::modulus(B) << endl;
        C=C-nint(scalar_product(B,C)/aurostd::modulus(B)/aurostd::modulus(C))*B;
        if(aurostd::modulus(C) > aurostd::modulus(B)-eps) {
          if(aurostd::mod(SwapCnt,2)!=0) {temp=B;B=C;C=temp;} // GH CORRECT
          // BUG BUG BUG if(aurostd::mod(it,2)!=0) {temp=B;B=C;C=temp;} // GH
          // GH Make sure the routine doesn't change the order of B and C on output
          // In other words, switch B and C again if odd number of swaps so far
          return TRUE; // basis cannot be further reduced
        }
      }
    }


  // ***************************************************************************************************
  // This routine takes three vectors, A,B,C, defining a lattice, and reduces the first one so that it
  // is a close as possible to the origin while remaining in the affine plane which is defined by B,C but
  // shifted by A. See Lecture notes in computer science, ISSN 0302-974, ANTS - VI : algorithmic
  // number theory, 2004, vol. 3076, pp. 338-357 ISBN 3-540-22156-5
  template<class utype> void
    reduce_A_in_ABC(xvector<utype>& A, xvector<utype>& B, xvector<utype>& C,utype eps) {
      xvector<utype> T(3);  // closest point to origin in B,C+A affine plane
      xmatrix<utype> ABC(3,3),ABCinv(3,3),oldABC(3,3); // Matrices of ABC basis vectors and inverse
      xvector<utype> dist(4); // the distances from T to enclosing lattice points of B,C (4 corners of the ppiped)
      xvector<utype> i(3),i1(3),i2(3),i3(3),i4(3); // lattice coordinates of A, in the affine plane, using the B,C basis vectors
      int idx; // index of the smallest distance from T to a lattice point in B,C
      //[CO20191201 - OBSOLETE]bool err;
      //integer j
      utype lambda;
      //print *,"entering reduction routine..."
      //write(*,'("aurostd::modulus(A): ",f7.3,5x," A ",3(f7.3,1x)A)') aurostd::modulus(A), A
      for(int i=1;i<=3;i++) {
        ABC(i,1)=A(i);
        ABC(i,2)=B(i);
        ABC(i,3)=C(i);
      }
      oldABC=ABC;
      // Use Gaussian reduction to reduce the B,C 2D basis so that it is itself Minkowski reduced. If this
      // is done then the closest lattice point (in B,C plane) to projection of A (into the B,C plane) is
      // guaranteed to be one of the corners of the unit cell enclosing the projection of A
      gaussian_reduce_two_vectors(B,C,eps);

      //do j=1,3
      //   write(*,'(3(f11.5,1x))') ABC(j,:)
      //enddo
      //
      // First thing to do is find the (real, not lattice) point in the affine plane B,C + A that is
      // nearest the origin. Call this T.
      lambda=-scalar_product(A,vector_product(B,C))/(aurostd::modulus(vector_product(B,C))*aurostd::modulus(vector_product(B,C)));
      T=A + lambda*vector_product(B,C);

      //print *,lambda
      //write(*,'("T (vec in B,C affine plane): ",3(f10.3,1x))') T

      // Now find the four points of the B,C lattice, in the affine plane, that enclose the point T
      for(int i=1;i<=3;i++) {//GH We need these 3 lines to load matrix ABC again with the vectors A,B,C
        ABC(i,1)=A(i);ABC(i,2)=B(i);ABC(i,3)=C(i);//GH
      }//GH
      if(aurostd::isNonInvertible(ABC)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"A,B,C vectors in reduce_A_in_ABC are co-planar",_VALUE_RANGE_);} //CO20191201
      ABCinv=inverse(ABC);  //CO20191201
      i=aurostd::nint(ABCinv*T);

      // print *,"Lattice coordinates of origin enclosing T:", i
      // Compute the distance from T to each of the four points and pick the one that is the closest.
      i1(1)=i(1);i1(2)=i(2);i1(3)=i(3);dist(1)=aurostd::modulus(T-ABC*i1);
      i2(1)=i(1);i2(2)=i(2)+1;i2(3)=i(3);dist(2)=aurostd::modulus(T-ABC*i2);
      i3(1)=i(1);i3(2)=i(2);i3(3)=i(3)+1;dist(3)=aurostd::modulus(T-ABC*i3);
      i4(1)=i(1);i4(2)=i(2)+1;i4(3)=i(3)+1;dist(4)=aurostd::modulus(T-ABC*i4);
      idx=0;idx=aurostd::mini(dist);
      //write(*,'("Dists: ",4(f10.5,1x))') dist

      //if(.not. equal(,origdist,eps)) then // Only change A if the new one really

      if(idx==1) A=A-ABC*i1;
      if(idx==2) A=A-ABC*i2;
      if(idx==3) A=A-ABC*i3;
      if(idx==4) A=A-ABC*i4;
      if(idx==0) {
        string message = "Case failed in reduce_A_in_ABC";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_RANGE_);
      }
      //endif
      //write(*,'("aurostd::modulus(A): ",f7.3,5x," A ",3(f7.3,1x)A)') aurostd::modulus(A), A
      for(int i=1;i<=3;i++) {//GH We need these 3 lines to load matrix ABC again with the vectors A,B,C
        ABC(i,1)=A(i);ABC(i,2)=B(i);ABC(i,3)=C(i);//GH
      }//GH

      // [OBSOLETE]  aurostd::matrix_inverse(ABC,ABCinv,err);
      //[CO20191201 - OBSOLETE]err=aurostd::inverse(ABC,ABCinv);
      //
      if(aurostd::isNonInvertible(ABC)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"A,B,C vectors in reduce_A_in_ABC are co-planar",_VALUE_RANGE_);} //CO20191201
      ABCinv=inverse(ABC);  //CO20191201
      ABC=ABCinv*oldABC-aurostd::nint(ABCinv*oldABC);

      for(int i=1;i<=3;i++) {
        for(int j=1;j<=3;j++) {
          if(aurostd::abs(ABC(i,j))>eps) {
            stringstream message;
            message << "eps=" << eps << std::endl;
            message << "ABC(i,j)=" << ABC(i,j) << std::endl;
            message << "ABCinv=" << ABCinv << std::endl;
            message << "oldABC=" << oldABC << std::endl;
            message << "ABCinv*oldABC=" << ABCinv*oldABC << std::endl;
            message << "ABCinv*oldABC-aurostd::nint(ABCinv*oldABC)=" << ABCinv*oldABC-aurostd::nint(ABCinv*oldABC) << std::endl;
            message << "Lattice was not preserved  in reduce_A_in_ABC";
            throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
          }
        }
      }
      //read(*,*)
    }

  // ***************************************************************************************************
  //  This routine takes a set of basis vectors (that form a lattice) and reduces them so that they form
  //  the shortest possible basis. See Lecture notes in computer science, ISSN 0302-974, ANTS - VI : algorithmic
  //  number theory, 2004, vol. 3076, pp. 338-357 ISBN 3-540-22156-5
  template<class utype> utype
    reduce_to_shortest_basis(const xmatrix<utype>& IN,xmatrix<utype>& OUT,utype eps,bool VERBOSE) {
      xvector<utype> A(3),B(3),C(3);
      xmatrix<utype> check(3,3);
      //[CO20191201 - OBSOLETE]bool err;
      utype od,odnew;
      int ii=0,iimax=10000;
      // IN has colum-vectors
      for(int i=1;i<=3;i++) {
        A(i)=IN(i,1); // 1st vector
        B(i)=IN(i,2); // 2nd vector
        C(i)=IN(i,3); // 3rd vector
      }
      odnew=orthogonality_defect(IN);
      if(VERBOSE) cout << "aurostd::reduce_to_shortest_basis: Before reduction, the orthogonality defect of the basis was " << odnew << endl;
      bool goexit=FALSE;
      while(goexit==FALSE) {
        od=odnew;
        reduce_A_in_ABC(A,B,C,eps);
        reduce_A_in_ABC(B,C,A,eps);
        reduce_A_in_ABC(C,A,B,eps);
        for(int i=1;i<=3;i++) {
          OUT(i,1)=A(i);OUT(i,2)=B(i);OUT(i,3)=C(i);
        }
        odnew=orthogonality_defect(OUT);
        // write(*,'("OD: ",2(f7.3,1x))') odnew, od
        //      cerr << od << " " << odnew << endl;
        if(aurostd::abs(od-odnew)<eps) goexit=TRUE;
        if(ii++>iimax) goexit=TRUE;
        //     if(ii++>iimax) {OUT=IN;return orthogonality_defect(OUT);}
      }
      if(aurostd::isNonInvertible(OUT)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"OUT matrix is singular in reduce_to_shortest_basis",_VALUE_RANGE_);} //CO20191201
      check=aurostd::inverse(OUT);
      //  Check that the conversion from old to new lattice vectors is still an integer matrix
      if(sum(abs(check*IN-nint(check*IN)))>eps) {
        string message = "Reduced lattice vectors in reduce_to_shortest_basis changed the original lattice";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(VERBOSE) cout << "aurostd::reduce_to_shortest_basis: After reduction, the orthogonality defect of the basis is " << orthogonality_defect(OUT) << endl;
      //GH if we have a left-handed basis, then exchange two vectors so that the basis is right-handed (I don't care but VASP does...Why?)
      if(aurostd::det(OUT)<eps) {
        utype temp;
        for(int i=1;i<=3;i++) {temp=OUT(i,1);OUT(i,1)=OUT(i,2);OUT(i,2)=temp;} // swap 1st with 2nd vector
      }
      // OUT has colum-vectors
      return orthogonality_defect(OUT);
    }

  template<class utype> xmatrix<utype>
    reduce_to_shortest_basis(const xmatrix<utype>& IN,utype eps,bool VERBOSE) {
      xmatrix<utype> newbasis(3,3);
      reduce_to_shortest_basis(IN,newbasis,eps,VERBOSE);
      return newbasis;
    }

  template<class utype> xmatrix<utype>
    reduce_to_shortest_basis(const xmatrix<utype>& IN) {
      return reduce_to_shortest_basis(IN,(utype)_DEFAULT_EPS_BASIS_REDUCE_,FALSE);
    }
}

//*****************************************************************************
// EIGENVECTORS EIGENVALUES STUFF

// ****************************************************
namespace aurostd {

#define NRANSI
#define ROTATE(a,i,j,k,l)   {g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);}
  template<class utype>
    int jacobi(const xmatrix<utype> &ain,xvector<utype> &d,xmatrix<utype> &v) {
      // Computes all eigenvalues and eigenvectors of a real symmetric xmatrix a[1..n][1..n].
      // On output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
      // v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
      // a. The function returns the number of Jacobi rotations that were required.

      stringstream message;
      int j,iq,ip,i,n,nrot=0;
      utype tresh,theta,tau,t,sm,s,h,g,c;
      xmatrix<utype> a(ain);
      n=a.rows;
      if(a.rows!=a.cols) {
        message << "'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(v.rows!=v.cols) {
        message << "'v' matrix not square  v.rows" << v.rows << " v.cols=" << v.cols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(a.rows!=v.rows) {
        message << "'a' and 'v' matrices must have same size  a.rows" << a.rows << " v.rows=" << v.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(a.rows!=d.rows) {
        message << "'a' and 'd' objects must have same size  a.rows" << a.rows << " d.rows=" << d.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      xvector<utype> b(1,n);
      xvector<utype> z(1,n);
      for (ip=1;ip<=n;ip++) {
        for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
        v[ip][ip]=1.0;
      }
      for (ip=1;ip<=n;ip++) {
        b[ip]=d[ip]=a[ip][ip];
        z[ip]=0.0;
      }
      nrot=0;
      for (i=1;i<=50;i++) {
        sm=0.0;
        for (ip=1;ip<=n-1;ip++) {
          for (iq=ip+1;iq<=n;iq++)
            sm += aurostd::abs(a[ip][iq]);
        }
        if(sm == 0.0) {
          //~z;~b;
          return nrot;
        }
        if(i < 4)
          tresh=0.2*sm/(n*n);
        else
          tresh=0.0;
        for (ip=1;ip<=n-1;ip++) {
          for (iq=ip+1;iq<=n;iq++) {
            g=100.0*aurostd::abs(a[ip][iq]);
            if(i > 4 && (utype)(aurostd::abs(d[ip])+g) == (utype)aurostd::abs(d[ip])
                && (utype)(aurostd::abs(d[iq])+g) == (utype)aurostd::abs(d[iq]))
              a[ip][iq]=0.0;
            else if(aurostd::abs(a[ip][iq]) > tresh) {
              h=d[iq]-d[ip];
              if((utype)(aurostd::abs(h)+g) == (utype)aurostd::abs(h)) {
                t=(a[ip][iq])/h;
              } else {
                theta=0.5*h/(a[ip][iq]);
                t=1.0/(aurostd::abs(theta)+aurostd::sqrt(1.0+theta*theta));
                if(theta < 0.0) t = -t;
              }
              c=1.0/sqrt(1+t*t);
              s=t*c;
              tau=s/(1.0+c);
              h=t*a[ip][iq];
              z[ip] -= h;
              z[iq] += h;
              d[ip] -= h;
              d[iq] += h;
              a[ip][iq]=0.0;
              for (j=1;j<=ip-1;j++) {ROTATE(a,j,ip,j,iq);}
              for (j=ip+1;j<=iq-1;j++) {ROTATE(a,ip,j,j,iq);}
              for (j=iq+1;j<=n;j++) {ROTATE(a,ip,j,iq,j);}
              for (j=1;j<=n;j++) {ROTATE(v,j,ip,j,iq);}
              ++(nrot);
            }
          }
        }
        for (ip=1;ip<=n;ip++) {
          b[ip] += z[ip];
          d[ip]=b[ip];
          z[ip]=0.0;
        }
      }
      throw aurostd::xerror(__AFLOW_FILE__, "xmatrix::jacobi()", "Too many iterations.", _RUNTIME_ERROR_);
    }
#undef ROTATE

  //ME20190815
  // Jacobi algorithm for Hermitian matrices (used to be in APL/apl_aplmath.cpp)
  // Based on http://arxiv.org/abs/physics/0607103
  template<class utype>
    xvector<utype> jacobiHermitian(xmatrix<xcomplex<utype> >& A, char _sort_) {
      xmatrix<xcomplex<utype> > U(A.rows, A.cols);
      return jacobiHermitian(A, U, _sort_);
    }

  template<class utype>
    xvector<utype> jacobiHermitian(xmatrix<xcomplex<utype> >& A, xmatrix<xcomplex<utype> >& U, char _sort_) {
      // Matrices have to be square
      if (!A.issquare) {
        string message = "Input matrix is not square.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      // Reshape eigenvector matrix if needed
      if ((U.rows != A.rows) || (U.cols != A.cols)) {
        xmatrix<xcomplex<utype> > U_new(A.urows, A.ucols, A.lrows, A.lcols);
        U = U_new;
      }

      // Initialize
      utype eps = std::numeric_limits<utype>::epsilon();
      uint max_sweeps = 50;
      uint threshold_sweep = 4;
      uint n = A.rows;
      double reduction = 0.04/std::pow(n, 4);

      xvector<utype> d(n);
      xmatrix<utype> ev(2, n);
      U.clear();
      // Initialize U as unit matrix and d as diagonal elements of A
      for (uint p = 1; p <= n; p++) {
        ev[1][p] = 0.0;
        ev[2][p] = A[p][p].re;
        d[p] = ev[2][p];
        U[p][p] = 1.0;
      }

      utype sum, threshold, t, delta, invc, s;
      xcomplex<utype> x, y, Apq, cApq;
      uint p, q, j, nSweep;

      // Perform sweeps
      for (nSweep = 1; nSweep <= max_sweeps; nSweep++) {
        // Convergence criterion: sum of the squares of the off-diagonal elements
        sum = 0.0;
        for (q = 2; q <= n; q++) {
          for (p = 1; p <= q - 1; p++) {
            sum += magsqr(A[p][q]);
          }
        }
        if (sum < 0.5 * eps) break;
        if (nSweep < threshold_sweep) threshold = reduction * sum;
        else threshold = 0;

        // Perform Jacobi rotations
        for (q = 2; q <= n; q++) {
          for (p = 1; p <= q - 1; p++) {
            sum = magsqr(A[p][q]);
            if ((nSweep > threshold_sweep) &&
                (sum < 0.5 * eps * std::max<utype>(ev[2][p] * ev[2][p], ev[2][q] * ev[2][q]))) {
              A[p][q] = 0;
            } else if (sum > threshold) {
              t = 0.5 * (ev[2][p] - ev[2][q]);
              t = 1.0/(t + copysign(sqrt(t * t + sum), t));
              delta = t * sum;
              ev[1][p] = ev[1][p] + delta;
              ev[2][p] = d[p] + ev[1][p];
              ev[1][q] = ev[1][q] - delta;
              ev[2][q] = d[q] + ev[1][q];

              invc = sqrt(delta * t + 1);
              s = t/invc;
              t = delta/(invc + 1);

              Apq = A[p][q];
              cApq = conj(Apq);

              for (j = 1; j <= p - 1; j++) {
                x = A[j][p];
                y = A[j][q];
                A[j][p] = x + s * (cApq * y - t * x);
                A[j][q] = y - s * (Apq * x + t * y);
              }

              for (j = p + 1; j <= q - 1; j++) {
                x = A[p][j];
                y = A[j][q];
                A[p][j] = x + s * (Apq * conj(y) - t * x);
                A[j][q] = y - s * (Apq * conj(x) + t * y);
              }

              for (j = q + 1; j <= n; j++) {
                x = A[p][j];
                y = A[q][j];
                A[p][j] = x + s * (Apq * y - t * x);
                A[q][j] = y - s * (cApq * x + t * y);
              }

              A[p][q] = 0;

              for (j = 1; j <= n; j++) {
                x = U[p][j];
                y = U[q][j];
                U[p][j] = x + s * (Apq * y - t * x);
                U[q][j] = y - s * (cApq * x + t * y);
              }
            }
          }
        }
        for (p = 1; p <= n; p++) {
          ev[1][p] = 0.0;
          d[p] = ev[2][p];
        }
      }

      if (nSweep > max_sweeps) {
        string message = "Number of sweeps exceeded maximum number of sweeps.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      // Sort - leave as is if sort mode not found
      if (_sort_ == 0) {  // ascending order
        utype temp;
        xcomplex<utype> xtemp;
        for (uint i = 1; i <= n - 1; i++) {
          for (uint j = i + 1; j <= n; j++) {
            if (d[j] < d[i]) {
              temp = d[j];
              d[j] = d[i];
              d[i] = temp;
              for (uint k = 1; k <= n; k++) {
                xtemp = U[j][k];
                U[j][k] = U[i][k];
                U[i][k] = xtemp;
              }
            }
          }
        }
      } else if (_sort_ == 1) {  // descending order
        utype temp;
        xcomplex<utype> xtemp;
        for (uint i = 1; i <= n - 1; i++) {
          for (uint j = i + 1; j <= n; j++) {
            if (d[j] > d[i]) {
              temp = d[j];
              d[j] = d[i];
              d[i] = temp;
              for (uint k = 1; k <= n; k++) {
                xtemp = U[j][k];
                U[j][k] = U[i][k];
                U[i][k] = xtemp;
              }
            }
          }
        }
      }

      // Transpose to have eigenvectors in columns
      U = trasp(U);
      return d;
    }
}

// ****************************************************
namespace aurostd {
  template<class utype>
    void eigsrt(xvector<utype> &d,xmatrix<utype> &v) {
      // Given the eigenvalues d[1..n]and eigenvectors v[1..n][1..n] as output fromjacobi
      // or tqli,this routine sorts the eigenvalues into descending order, and rearranges
      // the columns of v correspondingly. The method is straight insertion and is N2 rather than NlogN;
      // but since you have just done an N3 procedure to get the eigenvalues, you can afford yourself
      // this little indulgence.
      int k,j,i,n;
      utype p;

      n=v.rows;
      if(v.rows!=v.cols) {
        stringstream message;
        message << "'v' matrix not square  v.rows" << v.rows << " v.cols=" << v.cols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(v.rows!=d.rows) {
        stringstream message;
        message << "'v' and 'd' objects must have same size  v.rows" << v.rows << " d.rows=" << d.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      for (i=1;i<n;i++) {
        p=d[k=i];
        for (j=i+1;j<=n;j++)
          if(d[j] >= p) p=d[k=j];
        if(k != i) {
          d[k]=d[i];
          d[i]=p;
          for (j=1;j<=n;j++) {
            p=v[j][i];
            v[j][i]=v[j][k];
            v[j][k]=p;
          }
        }
      }
    }
}

// ****************************************************
//CO20171129
namespace aurostd {
  template<class utype>
    void QRDecomposition_HouseHolder(const xmatrix<utype>& mat_orig,xmatrix<utype>& Q,xmatrix<utype>& R,utype tol) {  //CO20191110
      return QRDecomposition_HouseHolder_MW(mat_orig,Q,R,tol);
      //_MW() and _TB() show to have about the same run time, but _MW() can be slightly faster
      //could be because of extra LDEBUG bool checks
      //_MW() might introduce more error into R as x spans full column everytime, but the error falls below 1e-15
      //_TB() uses more memory storing all v's
      //_MW() is shorter and easier to follow
    }
  template<class utype>
    void QRDecomposition_HouseHolder_MW(const xmatrix<utype>& mat_orig,xmatrix<utype>& Q,xmatrix<utype>& R,utype tol) {  //CO20191110
      // mat is mxn, m>=n
      // inspired by https://www.mathworks.com/matlabcentral/answers/169648-qr-factorization-using-householder-transformations
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " mat_orig=" << endl;cerr << mat_orig << endl;}
      if(mat_orig.rows<mat_orig.cols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m<n, please flip the matrix",_VALUE_ERROR_);}

      R=mat_orig; //reset

      utype vModulus;
      Q=eye<utype>(R.urows,R.urows,R.lrows,R.lrows); //reset
      for(int k=R.lcols;k<=R.ucols;k++) {
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " step k=" << k << endl;}
        xmatrix<utype> x(R.urows,R.lcols,R.lrows,R.lcols);
        x.setmat(R.getxmat(k,R.urows,k,k),k,R.lcols);  //x(k:m,1)=R(k:m,k);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " x=" << endl;cerr << x << endl;}
        xmatrix<utype> v(x);  //+x first
        v[k][v.lcols]=x[k][x.lcols]+aurostd::modulus(x);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " v(unnormalized)=" << endl;cerr << v << endl;}
        vModulus=aurostd::modulus(v);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " ||v||=" << vModulus << endl;}
        if(!iszero(vModulus,tol)){  //prevents division by 0
          v/=vModulus;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " v(  normalized)=" << endl;cerr << v << endl;}
          xmatrix<utype> u=(utype)2.0*trasp(R)*v;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " u=" << endl;cerr << u << endl;}
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " R( pre)=" << endl;cerr << R << endl;}
          R-=v*trasp(u);  //product HR
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " R(post)=" << endl;cerr << R << endl;}
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " Q( pre)=" << endl;cerr << Q << endl;}
          Q-=(utype)2.0*Q*v*trasp(v); //product QR
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " Q(post)=" << endl;cerr << Q << endl;}
        }
      }

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " mat_orig=" << endl;cerr << mat_orig << endl;
        cerr << __AFLOW_FUNC__ << " Q=" << endl;cerr << Q << endl;
        cerr << __AFLOW_FUNC__ << " R=" << endl;cerr << R << endl;
      }

      if(!aurostd::isequal(mat_orig,Q*R,tol)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"QR decomposition failed (A!=Q*R)",_RUNTIME_ERROR_);}
      if(!aurostd::isequal(trasp(Q)*Q,eye<utype>(Q.urows,Q.ucols,Q.lrows,Q.lcols),tol)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"QR decomposition failed (Q not orthonormal)",_RUNTIME_ERROR_);}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
    }
  template<class utype>
    void QRDecomposition_HouseHolder_TB(const xmatrix<utype>& mat_orig,xmatrix<utype>& Q,xmatrix<utype>& R,utype tol) {  //CO20191110
      // mat is mxn, m>=n
      // See Numerical Linear Algebra, Trefethen and Bau, pg. 73
      // this function stores household rotations (v) to create Q at the end
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " mat_orig=" << endl;cerr << mat_orig << endl;}
      if(mat_orig.rows<mat_orig.cols){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m<n, please flip the matrix",_VALUE_ERROR_);}

      R=mat_orig; //reset

      int i=0,k=0;
      xmatrix<utype> x,A; //since x and A changes size with each loop, let getxmatInPlace() handle it internally
      utype vModulus = (utype)0;
      std::vector<xmatrix<utype> > V; //we need to save v's, Q is calculated afterwards and needs all v's present
      for(k=R.lcols;k<=R.ucols;k++) {
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " step k=" << k << endl;}
        R.getxmatInPlace(x,k,R.urows,k,k);  //build R([k:m],l)
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " x=" << endl;cerr << x << endl;}
        //v_k=sign(x1)||x||e1+x
        xmatrix<utype> v(x);  //+x first
        v[v.lrows][v.lcols]+=aurostd::sign(x[x.lrows][x.lcols])*aurostd::modulus(x);  //only applies to first entry of v (e1)
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " v(unnormalized)=" << endl;cerr << v << endl;}
        vModulus=aurostd::modulus(v);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " ||v||=" << vModulus << endl;}
        if(!iszero(vModulus,tol)){  //prevents division by 0
          v/=vModulus;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " v(  normalized)=" << endl;cerr << v << endl;}
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " R( pre)=" << endl;cerr << R << endl;}
          R.getxmatInPlace(A,k,R.urows,k,R.ucols);  //build R([k:m],[k:m])
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " A( pre)=" << endl;cerr << A << endl;}
          A-=(utype)2.0*v*trasp(v)*A;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " A(post)=" << endl;cerr << A << endl;}
          R.setmat(A,k,k);  //store A back into R
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " R(post)=" << endl;cerr << R << endl;}
        }
        V.push_back(v);
      }

      Q=xmatrix<utype>(R.urows,R.urows,R.lrows,R.lrows); //reset
      xmatrix<utype> ek(R.urows,R.lcols,R.lrows,R.lcols); //create identity matrix column vector
      for(k=R.lcols;k<=R.urows;k++){  //calculate Q*e1,Q*e2...
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " Q( pre)=" << endl;cerr << Q << endl;}
        for(i=R.lrows;i<=R.urows;i++){ek[i][ek.lcols]=(i==k) ? (utype)1 : (utype)0;}
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " ek( pre)=" << endl;cerr << ek << endl;}
        for(i=R.ucols;i>=R.lcols;i--){
          ek.getxmatInPlace(x,i,R.urows,R.lcols,R.lcols);
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " x( pre)=" << endl;cerr << x << endl;}
          x-=(utype)2.0*V[i-R.lcols]*trasp(V[i-R.lcols])*x;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " x(post)=" << endl;cerr << x << endl;}
          ek.setmat(x,i,R.lcols);
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " ek(post)=" << endl;cerr << ek << endl;}
        Q.setmat(ek,R.lrows,k);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " Q(post)=" << endl;cerr << Q << endl;}
      }

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " mat_orig=" << endl;cerr << mat_orig << endl;
        cerr << __AFLOW_FUNC__ << " Q=" << endl;cerr << Q << endl;
        cerr << __AFLOW_FUNC__ << " R=" << endl;cerr << R << endl;
      }

      if(!aurostd::isequal(mat_orig,Q*R,tol)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"QR decomposition failed (A!=Q*R)",_RUNTIME_ERROR_);}
      if(!aurostd::isequal(trasp(Q)*Q,eye<utype>(Q.urows,Q.ucols,Q.lrows,Q.lcols),tol)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"QR decomposition failed (Q not orthonormal)",_RUNTIME_ERROR_);}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
    }
  template<class utype>
    void getEHermite(utype a,utype b,xmatrix<utype>& ehermite){ //CO+YL20191201
      //implementation is inspired by that found here: http://pydoc.net/GBpy/0.1.1/GBpy.tools/
      //original license: GNU-GPL Style.
      //Elementary Hermite transformation.
      //For integers a and b, E = ehermite(a,b) returns
      //an integer matrix with determinant 1 such that E * [a;b] = [g;0],
      //where g is the gcd of a and b.
      //E = ehermite(a,b)
      //This function is in some ways analogous to GIVENS.

      bool LDEBUG=(FALSE || XHOST.DEBUG);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

      utype gcd=0,x=0,y=0;
      GCD(a,b,gcd,x,y);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " gcd(" << a << "," << b << ")=" << gcd << ", x=" << x << ", y=" << y << endl;}
      //ehermite is 2x2
      if(ehermite.rows!=2 || ehermite.cols!=2){xmatrix<utype> ehermite_tmp(2,2);ehermite=ehermite_tmp;}
      if(gcd){
        ehermite[ehermite.lrows][ehermite.lcols]=x;
        ehermite[ehermite.lrows][ehermite.ucols]=y;       //urows=lrows+1, ucols=lcols+1
        ehermite[ehermite.urows][ehermite.lcols]=-b/gcd;  //urows=lrows+1, ucols=lcols+1
        ehermite[ehermite.urows][ehermite.ucols]=a/gcd;   //urows=lrows+1, ucols=lcols+1
      }else{
        ehermite[ehermite.lrows][ehermite.lcols]=(utype)1;
        ehermite[ehermite.lrows][ehermite.ucols]=(utype)0;       //urows=lrows+1, ucols=lcols+1
        ehermite[ehermite.urows][ehermite.lcols]=(utype)0;       //urows=lrows+1, ucols=lcols+1
        ehermite[ehermite.urows][ehermite.ucols]=(utype)1;       //urows=lrows+1, ucols=lcols+1
      }

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
    }
  template<class utype>
    void getSmithNormalForm(const xmatrix<utype>& A_in,xmatrix<utype>& U_out,xmatrix<utype>& V_out,xmatrix<utype>& S_out,double tol){  //CO+YL20191201
      //implementation is inspired by that found here: http://pydoc.net/GBpy/0.1.1/GBpy.tools/
      //original license: GNU-GPL Style.
      //Smith normal form of an integer matrix.
      //[U,S,V] = smith(A) returns integer matrices U, S, and V such that
      //S = U*A*V (rotated from A=U*S*V')
      //S is diagonal and nonnegative, S(i,i) divides S(i+1,i+1) for all i,
      //det U =+-1, and det V =+-1.
      //This function is in some ways analogous to SVD.
      //This looks much like an SVD algorithm that first bidiagonalizes
      //A by Givens rotations and then chases zeros, except for
      //the construction of the 2 by 2 elementary transformation.
      //we work with doubles inside, return int matrices later

      bool LDEBUG=(FALSE || XHOST.DEBUG);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " A=" << endl;cerr << A_in << endl;}

      xmatrix<double> S=aurostd::xmatrixutype2double(A_in);S.shift(1,1); //algorithm depends on lrows==lcols==1

      int m=S.rows,n=S.cols;
      int min_mn=std::min(m,n);
      xmatrix<double> U=eye<double>(m),V=eye<double>(n);

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " bidiagonalizing S with elementary Hermite transforms" << endl;}

      xmatrix<double> E(2,2);
      xmatrix<double> mXtwo_in(m,2),mXtwo_out(m,2),nXtwo_in(n,2),nXtwo_out(n,2),twoXn_in(2,n),twoXn_out(2,n);
      int j=0,i=0,jj=0;
      for(j=S.lcols;j<=min_mn;j++){
        //Zero column j below the diagonal.
        for(i=j+1;i<=m;i++){
          if(!iszero(S[i][j],tol)){
            //Construct an elementary Hermite transformation E
            //to zero S(i,j) by combining rows i and j.
            getEHermite(S[j][j],S[i][j],E);
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " getEHermite(S[j=" << j <<"][j=" << j << "]="<< S[j][j] <<",S[i=" << i << "][j=" << j << "]=" << S[i][j] << ")=" << endl;cerr << E << endl;}

            //Apply the transform to S
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(pre)=" << endl;cerr << S << endl;}
            //build S([j i],:)
            for(jj=twoXn_in.lcols;jj<=twoXn_in.ucols;jj++){
              twoXn_in[twoXn_in.lrows][jj]=S[j][jj];
              twoXn_in[twoXn_in.urows][jj]=S[i][jj];  //urows=lrows+1, ucols=lcols+1
            }
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " S([j=" << j << " i=" << i << "],:)=" << endl;cerr << twoXn_in << endl;}
            twoXn_out=E*twoXn_in; //2x2 x 2x3 = 2x3
            //store twoXn_out into S([j i],:)
            for(jj=twoXn_out.lcols;jj<=twoXn_out.ucols;jj++){
              S[j][jj]=twoXn_out[twoXn_in.lrows][jj];
              S[i][jj]=twoXn_out[twoXn_in.urows][jj];  //urows=lrows+1, ucols=lcols+1
            }
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(post)=" << endl;cerr << S << endl;}

            //Apply the transform to U
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(pre)=" << endl;cerr << U << endl;}
            //build U(:,[j i])
            for(jj=mXtwo_in.lrows;jj<=mXtwo_in.urows;jj++){
              mXtwo_in[jj][mXtwo_in.lcols]=U[jj][j];
              mXtwo_in[jj][mXtwo_in.ucols]=U[jj][i];  //urows=lrows+1, ucols=lcols+1
            }
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(:,[j=" << j << " i=" << i << "])=" << endl;cerr << mXtwo_in << endl;}
            mXtwo_out=mXtwo_in/E;
            //store mXtwo_out into U(:,[j i])
            for(jj=mXtwo_out.lrows;jj<=mXtwo_out.urows;jj++){
              U[jj][j]=mXtwo_out[jj][mXtwo_out.lcols];
              U[jj][i]=mXtwo_out[jj][mXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
            }
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(post)=" << endl;cerr << U << endl;}
          }
        }
        //Zero row j after the superdiagonal.
        for(i=j+2;i<=n;i++){
          if(!iszero(S[j][i],tol)){
            //Construct an elementary Hermite transformation E
            //to zero S(j,i) by combining columns j+1 and i.
            getEHermite(S[j][j+1],S[j][i],E);
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " getEHermite(S[j=" << j <<"][j+1=" << j+1 << "]="<< S[j][j+1] <<",S[j=" << j << "][i=" << i << "]=" << S[j][i] << ")=" << endl;cerr << E << endl;}

            //Apply the transform to S
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(pre)=" << endl;cerr << S << endl;}
            //build S(:,[j+1 i])
            for(jj=mXtwo_in.lrows;jj<=mXtwo_in.urows;jj++){
              mXtwo_in[jj][mXtwo_in.lcols]=S[jj][j+1];
              mXtwo_in[jj][mXtwo_in.ucols]=S[jj][i];  //urows=lrows+1, ucols=lcols+1
            }
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(:,[j+1=" << j+1 << " i=" << i << "])=" << endl;cerr << mXtwo_in << endl;}
            mXtwo_out=mXtwo_in*trasp(E);
            //store mXtwo_out into S(:,[j+1 i])
            for(jj=mXtwo_out.lrows;jj<=mXtwo_out.urows;jj++){
              S[jj][j+1]=mXtwo_out[jj][mXtwo_out.lcols];
              S[jj][i]=mXtwo_out[jj][mXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
            }
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(post)=" << endl;cerr << S << endl;}

            //Apply the transform to V
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(pre)=" << endl;cerr << V << endl;}
            //build V(:,[j+1 i])
            for(jj=nXtwo_in.lrows;jj<=nXtwo_in.urows;jj++){
              nXtwo_in[jj][nXtwo_in.lcols]=V[jj][j+1];
              nXtwo_in[jj][nXtwo_in.ucols]=V[jj][i];  //urows=lrows+1, ucols=lcols+1
            }
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(:,[j+1=" << j+1 << " i=" << i << "])=" << endl;cerr << nXtwo_in << endl;}
            nXtwo_out=nXtwo_in/E;
            //store nXtwo_out into V(:,[j+1 i])
            for(jj=nXtwo_out.lrows;jj<=nXtwo_out.urows;jj++){
              V[jj][j+1]=nXtwo_out[jj][nXtwo_out.lcols];
              V[jj][i]=nXtwo_out[jj][nXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
            }
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(post)=" << endl;cerr << V << endl;}
          }
        }
      }

      //if results differ slightly from matlab, check _GCD() and enable matlab implementation for gcd(1,1)
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " U=" <<endl;cerr << U << endl;
        cerr << __AFLOW_FUNC__ << " V=" <<endl;cerr << V << endl;
        cerr << __AFLOW_FUNC__ << " S=" <<endl;cerr << S << endl;
      }

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " S is now upper bidiagonal, eliminating superdiagonal non-zeros" << endl;}

      xvector<double> D=S.getdiag(1);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " D=" << D << endl;}

      int k=0;
      double q=0.0;
      while(!iszero(D,tol)){
        //Start chasing bulge at first nonzero superdiagonal element.
        k=-1;
        for(i=D.lrows;i<=D.urows;i++){
          if(!iszero(D[i],tol)){k=i;break;}
        }
        //be careful, k refers to index of D, not S

        //To guarantee reduction in S(k,k), first make S(k,k) positive
        //and make S(k,k+1) nonnegative and less than S(k,k).
        if(std::signbit(S[k][k])){
          for(i=S.lcols;i<=S.ucols;i++){S[k][i]=-S[k][i];}
          for(i=U.lrows;i<=U.urows;i++){U[i][k]=-U[i][k];}
        }
        q=std::floor(S[k][k+1]/S[k][k]);
        E[1][1]=1;E[1][2]=0;
        E[2][1]=-q;E[2][2]=1;

        //Apply the transform to S
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(pre)=" << endl;cerr << S << endl;}
        //build S(:,[k k+1])
        for(jj=mXtwo_in.lrows;jj<=mXtwo_in.urows;jj++){
          mXtwo_in[jj][mXtwo_in.lcols]=S[jj][k];
          mXtwo_in[jj][mXtwo_in.ucols]=S[jj][k+1];  //urows=lrows+1, ucols=lcols+1
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(:,[k=" << k << " k+1=" << k+1 << "])=" << endl;cerr << mXtwo_in << endl;}
        mXtwo_out=mXtwo_in*trasp(E);
        //store mXtwo_out into S(:,[k k+1])
        for(jj=mXtwo_out.lrows;jj<=mXtwo_out.urows;jj++){
          S[jj][k]=mXtwo_out[jj][mXtwo_out.lcols];
          S[jj][k+1]=mXtwo_out[jj][mXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(post)=" << endl;cerr << S << endl;}

        //Apply the transform to V
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(pre)=" << endl;cerr << V << endl;}
        //build V(:,[k k+1])
        for(jj=nXtwo_in.lrows;jj<=nXtwo_in.urows;jj++){
          nXtwo_in[jj][nXtwo_in.lcols]=V[jj][k];
          nXtwo_in[jj][nXtwo_in.ucols]=V[jj][k+1];  //urows=lrows+1, ucols=lcols+1
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(:,[k=" << k << " k+1=" << k+1 << "])=" << endl;cerr << nXtwo_in << endl;}
        nXtwo_out=nXtwo_in/E;
        //store nXtwo_out into V(:,[k k+1])
        for(jj=nXtwo_out.lrows;jj<=nXtwo_out.urows;jj++){
          V[jj][k]=nXtwo_out[jj][nXtwo_out.lcols];
          V[jj][k+1]=nXtwo_out[jj][nXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(post)=" << endl;cerr << V << endl;}

        if(!iszero(S[k][k+1],tol)){
          //Zero the first nonzero superdiagonal element
          //using columns k and k+1, to start the bulge at S(k+1,k).
          getEHermite(S[k][k],S[k][k+1],E);
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " getEHermite(S[k=" << k <<"][k=" << k << "]="<< S[k][k] <<",S[k=" << k << "][k+1=" << k+1 << "]=" << S[k][k+1] << ")=" << endl;cerr << E << endl;}

          //Apply the transform to S
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(pre)=" << endl;cerr << S << endl;}
          //build S(:,[k k+1])
          for(jj=mXtwo_in.lrows;jj<=mXtwo_in.urows;jj++){
            mXtwo_in[jj][mXtwo_in.lcols]=S[jj][k];
            mXtwo_in[jj][mXtwo_in.ucols]=S[jj][k+1];  //urows=lrows+1, ucols=lcols+1
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(:,[k=" << k << " k+1=" << k+1 << "])=" << endl;cerr << mXtwo_in << endl;}
          mXtwo_out=mXtwo_in*trasp(E);
          //store mXtwo_out into S(:,[k k+1])
          for(jj=mXtwo_out.lrows;jj<=mXtwo_out.urows;jj++){
            S[jj][k]=mXtwo_out[jj][mXtwo_out.lcols];
            S[jj][k+1]=mXtwo_out[jj][mXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(post)=" << endl;cerr << S << endl;}

          //Apply the transform to V
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(pre)=" << endl;cerr << V << endl;}
          //build V(:,[k k+1])
          for(jj=nXtwo_in.lrows;jj<=nXtwo_in.urows;jj++){
            nXtwo_in[jj][nXtwo_in.lcols]=V[jj][k];
            nXtwo_in[jj][nXtwo_in.ucols]=V[jj][k+1];  //urows=lrows+1, ucols=lcols+1
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(:,[k=" << k << " k+1=" << k+1 << "])=" << endl;cerr << nXtwo_in << endl;}
          nXtwo_out=nXtwo_in/E;
          //store nXtwo_out into V(:,[k k+1])
          for(jj=nXtwo_out.lrows;jj<=nXtwo_out.urows;jj++){
            V[jj][k]=nXtwo_out[jj][nXtwo_out.lcols];
            V[jj][k+1]=nXtwo_out[jj][nXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(post)=" << endl;cerr << V << endl;}

          for(j=S.lcols;j<=min_mn;j++){
            if(j+1<=m){
              //Zero S(j+1,j) using rows j and j+1.
              getEHermite(S[j][j],S[j+1][j],E);
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " getEHermite(S[j=" << j <<"][j=" << j << "]="<< S[j][j] <<",S[j+1=" << j+1 << "][j=" << j << "]=" << S[j+1][j] << ")=" << endl;cerr << E << endl;}

              //Apply the transform to S
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(pre)=" << endl;cerr << S << endl;}
              //build S([j j+1],:)
              for(jj=twoXn_in.lcols;jj<=twoXn_in.ucols;jj++){
                twoXn_in[twoXn_in.lrows][jj]=S[j][jj];
                twoXn_in[twoXn_in.urows][jj]=S[j+1][jj];  //urows=lrows+1, ucols=lcols+1
              }
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " S([j=" << j << " j+1=" << j+1 << "],:)=" << endl;cerr << twoXn_in << endl;}
              twoXn_out=E*twoXn_in; //2x2 x 2x3 = 2x3
              //store twoXn_out into S([j j+1],:)
              for(jj=twoXn_out.lcols;jj<=twoXn_out.ucols;jj++){
                S[j][jj]=twoXn_out[twoXn_in.lrows][jj];
                S[j+1][jj]=twoXn_out[twoXn_in.urows][jj];  //urows=lrows+1, ucols=lcols+1
              }
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(post)=" << endl;cerr << S << endl;}

              //Apply the transform to U
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(pre)=" << endl;cerr << U << endl;}
              //build U(:,[j j+1])
              for(jj=mXtwo_in.lrows;jj<=mXtwo_in.urows;jj++){
                mXtwo_in[jj][mXtwo_in.lcols]=U[jj][j];
                mXtwo_in[jj][mXtwo_in.ucols]=U[jj][j+1];  //urows=lrows+1, ucols=lcols+1
              }
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(:,[j=" << j << " j+1=" << j+1 << "])=" << endl;cerr << mXtwo_in << endl;}
              mXtwo_out=mXtwo_in/E;
              //store mXtwo_out into U(:,[j j+1])
              for(jj=mXtwo_out.lrows;jj<=mXtwo_out.urows;jj++){
                U[jj][j]=mXtwo_out[jj][mXtwo_out.lcols];
                U[jj][j+1]=mXtwo_out[jj][mXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
              }
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(post)=" << endl;cerr << U << endl;}
            }
            if(j+2<=n){
              //Zero S(j,j+2) using columns j+1 and j+2.
              getEHermite(S[j][j+1],S[j][j+2],E);
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " getEHermite(S[j=" << j <<"][j+1=" << j+1 << "]="<< S[j][j+1] <<",S[j=" << j << "][j+2=" << j+2 << "]=" << S[j][j+2] << ")=" << endl;cerr << E << endl;}

              //Apply the transform to S
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(pre)=" << endl;cerr << S << endl;}
              //build S(:,[j+1 j+2])
              for(jj=mXtwo_in.lrows;jj<=mXtwo_in.urows;jj++){
                mXtwo_in[jj][mXtwo_in.lcols]=S[jj][j+1];
                mXtwo_in[jj][mXtwo_in.ucols]=S[jj][j+2];  //urows=lrows+1, ucols=lcols+1
              }
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(:,[j+1=" << j+1 << " j+2=" << j+2 << "])=" << endl;cerr << mXtwo_in << endl;}
              mXtwo_out=mXtwo_in*trasp(E);
              //store mXtwo_out into S(:,[j+1 j+2])
              for(jj=mXtwo_out.lrows;jj<=mXtwo_out.urows;jj++){
                S[jj][j+1]=mXtwo_out[jj][mXtwo_out.lcols];
                S[jj][j+2]=mXtwo_out[jj][mXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
              }
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(post)=" << endl;cerr << S << endl;}

              //Apply the transform to V
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(pre)=" << endl;cerr << V << endl;}
              //build V(:,[j+1 j+2])
              for(jj=nXtwo_in.lrows;jj<=nXtwo_in.urows;jj++){
                nXtwo_in[jj][nXtwo_in.lcols]=V[jj][j+1];
                nXtwo_in[jj][nXtwo_in.ucols]=V[jj][j+2];  //urows=lrows+1, ucols=lcols+1
              }
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(:,[j+1=" << j+1 << " j+2=" << j+2 << "])=" << endl;cerr << nXtwo_in << endl;}
              nXtwo_out=nXtwo_in/E;
              //store nXtwo_out into V(:,[j+1 j+2])
              for(jj=nXtwo_out.lrows;jj<=nXtwo_out.urows;jj++){
                V[jj][j+1]=nXtwo_out[jj][nXtwo_out.lcols];
                V[jj][j+2]=nXtwo_out[jj][nXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
              }
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(post)=" << endl;cerr << V << endl;}
            }
          }
        }
        D=S.getdiag(1);
      }

      //if results differ slightly from matlab, check _GCD() and enable matlab implementation for gcd(1,1)
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " U=" <<endl;cerr << U << endl;
        cerr << __AFLOW_FUNC__ << " V=" <<endl;cerr << V << endl;
        cerr << __AFLOW_FUNC__ << " S=" <<endl;cerr << S << endl;
      }

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " S is now diagonal, make it non-negative" << endl;}

      for(j=S.lcols;j<=min_mn;j++){
        if(std::signbit(S[j][j])){
          for(i=S.lcols;i<=S.ucols;i++){S[j][i]=-S[j][i];}
          for(i=U.lrows;i<=U.urows;i++){U[i][j]=-U[i][j];}
        }
      }

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " U=" <<endl;cerr << U << endl;
        cerr << __AFLOW_FUNC__ << " V=" <<endl;cerr << V << endl;
        cerr << __AFLOW_FUNC__ << " S=" <<endl;cerr << S << endl;
      }

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " squeezing factors to lower right to enforce divisibility condition" << endl;}

      double a=0.0,b=0.0,gcd=0.0,x=0.0,y=0.0;
      xmatrix<double> F(E),twoXtwo_in(2,2),twoXtwo_out(2,2);
      for(i=S.lcols;i<=min_mn;i++){
        for(j=i+1;j<=min_mn;j++){
          //Replace S(i,i), S(j,j) by their gcd and lcm respectively.
          a=S[i][i];
          b=S[j][j];
          GCD(a,b,gcd,x,y);

          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " a=" << a << endl;
            cerr << __AFLOW_FUNC__ << " b=" << b << endl;
            cerr << __AFLOW_FUNC__ << " gcd=" << gcd << endl;
            cerr << __AFLOW_FUNC__ << " x=" << x << endl;
            cerr << __AFLOW_FUNC__ << " y=" << y << endl;
          }

          E[1][1]=1.0;E[1][2]=y;
          E[2][1]=-b/gcd;E[2][2]=a*x/gcd;

          F[1][1]=x;F[1][2]=1.0;
          F[2][1]=-b*y/gcd;F[2][2]=a/gcd;

          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " E=" << endl;cerr << E << endl;
            cerr << __AFLOW_FUNC__ << " F=" << endl;cerr << F << endl;
          }

          //Apply the transform to S
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(pre)=" << endl;cerr << S << endl;}
          //build S([i j],[i j])
          twoXtwo_in[twoXtwo_in.lrows][twoXtwo_in.lcols]=S[i][i];twoXtwo_in[twoXtwo_in.lrows][twoXtwo_in.ucols]=S[i][j];  //urows=lrows+1, ucols=lcols+1
          twoXtwo_in[twoXtwo_in.urows][twoXtwo_in.lcols]=S[j][i];twoXtwo_in[twoXtwo_in.urows][twoXtwo_in.ucols]=S[j][j];  //urows=lrows+1, ucols=lcols+1
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " S([i=" << i << " j=" << j << "],[i=" << i << " j=" << j << "])=" << endl;cerr << twoXtwo_in << endl;}
          twoXtwo_out=E*twoXtwo_in*trasp(F);
          //store twoXtwo_out into S([i j],[i j])
          S[i][i]=twoXtwo_out[twoXtwo_out.lrows][twoXtwo_out.lcols];S[i][j]=twoXtwo_out[twoXtwo_out.lrows][twoXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
          S[j][i]=twoXtwo_out[twoXtwo_out.urows][twoXtwo_out.lcols];S[j][j]=twoXtwo_out[twoXtwo_out.urows][twoXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " S(post)=" << endl;cerr << S << endl;}

          //Apply the transform to U
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(pre)=" << endl;cerr << U << endl;}
          //build U(:,[i j])
          for(jj=mXtwo_in.lrows;jj<=mXtwo_in.urows;jj++){
            mXtwo_in[jj][mXtwo_in.lcols]=U[jj][i];
            mXtwo_in[jj][mXtwo_in.ucols]=U[jj][j];  //urows=lrows+1, ucols=lcols+1
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(:,[i=" << i << " j=" << j << "])=" << endl;cerr << mXtwo_in << endl;}
          mXtwo_out=mXtwo_in/E;
          //store mXtwo_out into U(:,[i j])
          for(jj=mXtwo_out.lrows;jj<=mXtwo_out.urows;jj++){
            U[jj][i]=mXtwo_out[jj][mXtwo_out.lcols];
            U[jj][j]=mXtwo_out[jj][mXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " U(post)=" << endl;cerr << U << endl;}

          //Apply the transform to V
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(pre)=" << endl;cerr << V << endl;}
          //build V(:,[i j])
          for(jj=nXtwo_in.lrows;jj<=nXtwo_in.urows;jj++){
            nXtwo_in[jj][nXtwo_in.lcols]=V[jj][i];
            nXtwo_in[jj][nXtwo_in.ucols]=V[jj][j];  //urows=lrows+1, ucols=lcols+1
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(:,[i=" << i << " j=" << j << "])=" << endl;cerr << nXtwo_in << endl;}
          nXtwo_out=nXtwo_in/F;
          //store nXtwo_out into V(:,[i j])
          for(jj=nXtwo_in.lrows;jj<=nXtwo_in.urows;jj++){
            V[jj][i]=nXtwo_out[jj][nXtwo_out.lcols];
            V[jj][j]=nXtwo_out[jj][nXtwo_out.ucols];  //urows=lrows+1, ucols=lcols+1
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " V(post)=" << endl;cerr << V << endl;}
        }
      }

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " U=" <<endl;cerr << U << endl;
        cerr << __AFLOW_FUNC__ << " V=" <<endl;cerr << V << endl;
        cerr << __AFLOW_FUNC__ << " S=" <<endl;cerr << S << endl;
      }

      //CONVERT TO INTEGERS FIRST!
      //inverse of an integer matrix is an integer matrix IFF det(M)= 1/-1 (true for V and U as above)
      //algorithm is MUCH more stable this way
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " converting to xmatrix<int>" << endl;}
      U_out=xmatrixdouble2utype<utype>(U);
      V_out=xmatrixdouble2utype<utype>(V);
      S_out=xmatrixdouble2utype<utype>(S);
      //if results differ slightly from matlab, check _GCD() and enable matlab implementation for gcd(1,1)
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " U=" <<endl;cerr << U_out << endl;
        cerr << __AFLOW_FUNC__ << " V=" <<endl;cerr << V_out << endl;
        cerr << __AFLOW_FUNC__ << " S=" <<endl;cerr << S_out << endl;
      }

      //the routine should give SNF such that A=U*S*V'
      //we will rotate after this if the Matlab output is desired
      if(!aurostd::isequal(A_in,U_out*S_out*trasp(V_out),(utype)tol)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"SNF decomposition failed",_RUNTIME_ERROR_);}

      //operations below here for Matlab-like output

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " transposing V to match Matlab" << endl;}
      traspInPlace(V_out);  //such that A=U*S*V  and not A=U*S*V'
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " V=" <<endl;cerr << V_out << endl;}

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " inverting V and U to match Matlab" << endl;}
      V_out=inverse(V_out);U_out=inverse(U_out);  //to be identical to matlab's smithForm we need V -> inv(V) U-> inv(U)
      //if results differ slightly from matlab, check _GCD() and enable matlab implementation for gcd(1,1)
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " U=" <<endl;cerr << U_out << endl;
        cerr << __AFLOW_FUNC__ << " V=" <<endl;cerr << V_out << endl;
        cerr << __AFLOW_FUNC__ << " S=" <<endl;cerr << S_out << endl;
      }

      //shift everything to match A_in
      U_out.shift(A_in.lrows,A_in.lcols);
      V_out.shift(A_in.lrows,A_in.lcols);
      S_out.shift(A_in.lrows,A_in.lcols);

      //Matlab gives SNF such that S=U*A*V
      if(!aurostd::isequal(S_out,U_out*A_in*V_out,(utype)tol)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"SNF decomposition failed AFTER Matlab transformations",_RUNTIME_ERROR_);}

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
    }
}

// ****************************************************
namespace aurostd {
  template<class utype>
    void tred2(const xmatrix<utype> &a,xvector<utype> &d,xvector<utype> &e) {
      // Householder reduction of a real, symmetric matrix a[1..n][1..n].
      // On output, a is replaced by the orthogonal matrix Q eﬀecting the
      // transformation. d[1..n] returns the diagonal elments of
      // the tridiagonal matrix, and e[1..n] the oﬀ-diagonal elements, with e[1]=0.
      stringstream message;

      int l,k,j,i,n;
      utype scale,hh,h,g,f;

      n=a.rows;
      if(a.rows!=a.cols) {
        message << "'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(a.rows!=d.rows) {
        message << "'a' and 'd' objects must have same size  a.rows" << a.rows << " d.rows=" << d.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(a.rows!=e.rows) {
        message << "'a' and 'e' objects must have same size  a.rows" << a.rows << " e.rows=" << e.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      for (i=n;i>=2;i--) {
        l=i-1;
        h=scale=0.0;
        if(l > 1) {
          for (k=1;k<=l;k++)
            scale += aurostd::abs(a[i][k]);
          if(scale == 0.0)
            e[i]=a[i][l];
          else {
            for (k=1;k<=l;k++) {
              a[i][k] /= scale;
              h += a[i][k]*a[i][k];
            }
            f=a[i][l];
            g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i][l]=f-g;
            f=0.0;
            for (j=1;j<=l;j++) {
              a[j][i]=a[i][j]/h;
              g=0.0;
              for (k=1;k<=j;k++)
                g += a[j][k]*a[i][k];
              for (k=j+1;k<=l;k++)
                g += a[k][j]*a[i][k];
              e[j]=g/h;
              f += e[j]*a[i][j];
            }
            hh=f/(h+h);
            for (j=1;j<=l;j++) {
              f=a[i][j];
              e[j]=g=e[j]-hh*f;
              for (k=1;k<=j;k++)
                a[j][k] -= (f*e[k]+g*a[i][k]);
            }
          }
        } else
          e[i]=a[i][l];
        d[i]=h;
      }
      d[1]=0.0;
      e[1]=0.0;
      // Contents of this loop can be omitted if eigenvectors not
      // wanted except for statement d[i]=a[i][i];
      for (i=1;i<=n;i++) {
        l=i-1;
        if(d[i]) {
          for (j=1;j<=l;j++) {
            g=0.0;
            for (k=1;k<=l;k++)
              g += a[i][k]*a[k][j];
            for (k=1;k<=l;k++)
              a[k][j] -= g*a[k][i];
          }
        }
        d[i]=a[i][i];
        a[i][i]=1.0;
        for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
      }
    }

}

// ****************************************************
namespace aurostd {

  template<class utype>
    utype NR_SQR(utype a) {
      if(a==(utype) 0.0) return 0.0; else return a*a;
    }

  template<class utype>
    utype pythag(utype a, utype b) {
      utype absa,absb;
      absa=aurostd::abs(a);
      absb=aurostd::abs(b);
      if(absa > absb) return absa*sqrt(1.0+NR_SQR(absb/absa));
      else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+NR_SQR(absa/absb)));
    }


#define NR_SIGN(a,b) ((b) >= 0.0 ? aurostd::abs(a) : -aurostd::abs(a))
  template<class utype>
    void tqli(xvector<utype> &d,xvector<utype> &e,xmatrix<utype> &z) {
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
      int m,l,iter,i,k,n;
      utype s,r,p,g,f,dd,c,b;

      stringstream message;

      n=z.rows;
      if(z.rows!=z.cols) {
        message << "'z' matrix not square  z.rows" << z.rows << " z.cols=" << z.cols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(z.rows!=d.rows) {
        message << "'z' and 'd' objects must have same size  z.rows" << z.rows << " d.rows=" << d.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if(z.rows!=e.rows) {
        message << "'z' and 'e' objects must have same size  z.rows" << z.rows << " e.rows=" << e.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      for (i=2;i<=n;i++) e[i-1]=e[i];
      e[n]=0.0;
      for (l=1;l<=n;l++) {
        iter=0;
        do {
          for (m=l;m<=n-1;m++) {
            dd=aurostd::abs(d[m])+aurostd::abs(d[m+1]);
            if((utype)(aurostd::abs(e[m])+dd) == dd) break;
          }
          if(m != l) {
            if(iter++ == 30) {
              message << "Too many iterations in tqli.";
              throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
            }
            g=(d[l+1]-d[l])/(2.0*e[l]);
            r=pythag(g,(utype) 1.0);
            g=d[m]-d[l]+e[l]/(g+NR_SIGN(r,g));
            s=c=1.0;
            p=0.0;
            for (i=m-1;i>=l;i--) {
              f=s*e[i];
              b=c*e[i];
              e[i+1]=(r=pythag(f,g));
              if(r == 0.0) {
                d[i+1] -= p;
                e[m]=0.0;
                break;
              }
              s=f/r;
              c=g/r;
              g=d[i+1]-p;
              r=(d[i]-g)*s+2.0*c*b;
              d[i+1]=g+(p=s*r);
              g=c*r-b;
              for (k=1;k<=n;k++) {
                f=z[k][i+1];
                z[k][i+1]=s*z[k][i]+c*f;
                z[k][i]=c*z[k][i]-s*f;
              }
            }
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]=g;
            e[m]=0.0;
          }
        } while (m != l);
      }
    }
#undef NR_SIGN
}

// ****************************************************
namespace aurostd {
#define RADIX 2.0
  template<class utype>
    void balanc(xmatrix<utype> &a) {
      // Given a matrix a[1..n][1..n], this routine replaces it by
      // a balanced matrix with i dentical eigenvalues. A symmetric matrix
      // is already balanced and is unaﬀected by this procedure. The
      // parameter RADIX should be the machine’s ﬂoating-point radix.
      int last,j,i,n;
      utype s,r,g,f,c,sqrdx;

      n=a.rows;
      if(a.rows!=a.cols) {
        stringstream message;
        message << "'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      sqrdx=RADIX*RADIX;
      last=0;
      while (last == 0) {
        last=1;
        for (i=1;i<=n;i++) {
          r=c=0.0;
          for (j=1;j<=n;j++)
            if(j != i) {
              c += aurostd::abs(a[j][i]);
              r += aurostd::abs(a[i][j]);
            }
          if(c && r) {
            g=r/RADIX;
            f=1.0;
            s=c+r;
            while (c<g) {
              f *= RADIX;
              c *= sqrdx;
            }
            g=r*RADIX;
            while (c>g) {
              f /= RADIX;
              c /= sqrdx;
            }
            if((c+r)/f < 0.95*s) {
              last=0;
              g=1.0/f;
              for (j=1;j<=n;j++) a[i][j] *= g;
              for (j=1;j<=n;j++) a[j][i] *= f;
            }
          }
        }
      }
    }
}
#undef RADIX

// ****************************************************
namespace aurostd {
#define ELMHES_SWAP(g,h) {y=(g);(g)=(h);(h)=y;}
  template<class utype>
    // Reduction to Hessenberg form by the elimination method.
    // The real, nonsymmetric matrix a[1..n][1..n] is replaced by an upper
    // Hessenberg matrix with identical eigenvalues.
    // Recommended, but not required, is that this routine be preceded
    // by balanc. On output, the Hessenberg matrix is in elements a[i][j] with i<=j+1.
    // Elements with i>j+1 are to be thought of as zero,
    // but are returned with random values.
    void elmhes(xmatrix<utype> &a) {
      int m,j,i,n;
      utype y,x;

      n=a.rows;
      if(a.rows!=a.cols) {
        stringstream message;
        message << "'a' matrix not square  a.rows" << a.rows << " a.cols=" << a.cols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      for (m=2;m<n;m++) {
        x=0.0;
        i=m;
        for (j=m;j<=n;j++) {
          if(aurostd::abs(a[j][m-1]) > aurostd::abs(x)) {
            x=a[j][m-1];
            i=j;
          }
        }
        if(i != m) {
          for (j=m-1;j<=n;j++) ELMHES_SWAP(a[i][j],a[m][j]);
          for (j=1;j<=n;j++) ELMHES_SWAP(a[j][i],a[j][m]);
        }
        if(x) {
          for (i=m+1;i<=n;i++) {
            if((y=a[i][m-1]) != 0.0) {
              y /= x;
              a[i][m-1]=y;
              for (j=m;j<=n;j++)
                a[i][j] -= y*a[m][j];
              for (j=1;j<=n;j++)
                a[j][m] += y*a[j][i];
            }
          }
        }
      }
    }
#undef ELMHES_SWAP
}

// ****************************************************
namespace aurostd {
  // Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n].
  // On input a can be exactly as output from elmhes; on output it is destroyed.
  // The real and imaginary parts of the eigenvalues are returned in
  //wr[1..n] and wi[1..n], respectively.

#define NR_SIGN(a,b) ((b) >= 0.0 ? aurostd::abs(a) : -aurostd::abs(a))

  template<class utype>
    void hqr(xmatrix<utype> &a,xvector<utype> &wr,xvector<utype> &wi) {
      int nn,m,l,k,j,its,i,mmin,n;
      utype z,y,x,w,v,u,t,s,r=0,q=0,p=0,anorm;

      n=a.rows;
      if(a.rows!=a.cols) {
        stringstream message;
        message << "'a' matrix not square  a.rows=" << a.rows << " a.cols=" << a.cols;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      else if(wr.rows!=a.rows || wi.rows!=a.rows) { //SD20220420 - output lengths need to match input
        stringstream message;
        message << "'wr' and 'wi' lengths must be compatible with 'a' matrix, a.rows=" << a.rows << " wr.rows=" << wr.rows<< " wi.rows=" << wi.rows;
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      anorm=aurostd::abs(a[1][1]);
      for (i=2;i<=n;i++)
        for (j=(i-1);j<=n;j++)
          anorm += aurostd::abs(a[i][j]);
      nn=n;
      t=0.0;
      while (nn >= 1) {
        its=0;
        do {
          for (l=nn;l>=2;l--) {
            s=aurostd::abs(a[l-1][l-1])+aurostd::abs(a[l][l]);
            if(s == 0.0) s=anorm;
            if((utype)(aurostd::abs(a[l][l-1]) + s) == s) break;
          }
          x=a[nn][nn];
          if(l == nn) {
            wr[nn]=x+t;
            wi[nn--]=0.0;
          } else {
            y=a[nn-1][nn-1];
            w=a[nn][nn-1]*a[nn-1][nn];
            if(l == (nn-1)) {
              p=0.5*(y-x);
              q=p*p+w;
              z=sqrt(aurostd::abs(q));
              x += t;
              if(q >= 0.0) {
                z=p+NR_SIGN(z,p);
                wr[nn-1]=wr[nn]=x+z;
                if(z) wr[nn]=x-w/z;
                wi[nn-1]=wi[nn]=0.0;
              } else {
                wr[nn-1]=wr[nn]=x+p;
                wi[nn-1]= -(wi[nn]=z);
              }
              nn -= 2;
            } else {
              if(its == 30) {
                string message = "Too many iterations in hqr";
                throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
              }
              if(its == 10 || its == 20) {
                t += x;
                for (i=1;i<=nn;i++) a[i][i] -= x;
                s=aurostd::abs(a[nn][nn-1])+aurostd::abs(a[nn-1][nn-2]);
                y=x=0.75*s;
                w = -0.4375*s*s;
              }
              ++its;
              for (m=(nn-2);m>=l;m--) {
                z=a[m][m];
                r=x-z;
                s=y-z;
                p=(r*s-w)/a[m+1][m]+a[m][m+1];
                q=a[m+1][m+1]-z-r-s;
                r=a[m+2][m+1];
                s=aurostd::abs(p)+aurostd::abs(q)+aurostd::abs(r);
                p /= s;
                q /= s;
                r /= s;
                if(m == l) break;
                u=aurostd::abs(a[m][m-1])*(aurostd::abs(q)+aurostd::abs(r));
                v=aurostd::abs(p)*(aurostd::abs(a[m-1][m-1])
                    +aurostd::abs(z)+aurostd::abs(a[m+1][m+1]));
                if((utype)(u+v) == v) break;
              }
              for (i=m+2;i<=nn;i++) {
                a[i][i-2]=0.0;
                if(i != (m+2)) a[i][i-3]=0.0;
              }
              for (k=m;k<=nn-1;k++) {
                if(k != m) {
                  p=a[k][k-1];
                  q=a[k+1][k-1];
                  r=0.0;
                  if(k != (nn-1)) r=a[k+2][k-1];
                  if((x=aurostd::abs(p)+aurostd::abs(q)+aurostd::abs(r)) != 0.0) {
                    p /= x;
                    q /= x;
                    r /= x;
                  }
                }
                if((s=NR_SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
                  if(k == m) {
                    if(l != m)
                      a[k][k-1] = -a[k][k-1];
                  } else
                    a[k][k-1] = -s*x;
                  p += s;
                  x=p/s;
                  y=q/s;
                  z=r/s;
                  q /= p;
                  r /= p;
                  for (j=k;j<=nn;j++) {
                    p=a[k][j]+q*a[k+1][j];
                    if(k != (nn-1)) {
                      p += r*a[k+2][j];
                      a[k+2][j] -= p*z;
                    }
                    a[k+1][j] -= p*y;
                    a[k][j] -= p*x;
                  }
                  mmin = nn<k+3 ? nn : k+3;
                  for (i=l;i<=mmin;i++) {
                    p=x*a[i][k]+y*a[i][k+1];
                    if(k != (nn-1)) {
                      p += z*a[i][k+2];
                      a[i][k+2] -= p*r;
                    }
                    a[i][k+1] -= p*q;
                    a[i][k] -= p;
                  }
                }
              }
            }
          }
        } while (l < nn-1);
      }
    }
#undef NR_SIGN
}

// ****************************************************
namespace aurostd {
  // Finds all eigenvalues of matrix a[1..n][1..n]. The real and imaginary parts
  // of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.

  template<class utype>
    void eigen(const xmatrix<utype> &ain,xvector<utype> &wr,xvector<utype> &wi) {
      xmatrix<utype> a(ain);
      balanc(a);
      elmhes(a);
      hqr(a,wr,wi);
    }
}


//*****************************************************************************
// MATRIX NORMS - ME20190718
//*****************************************************************************

namespace aurostd {
  template<class utype>
    utype l1_norm(const xmatrix<utype>& m) {
      xvector<utype> vals(m.lcols, m.ucols);
      for (int i = m.lcols; i <= m.ucols; i++) {
        for (int j = m.lrows; j <= m.urows; j++) {
          vals[i] += abs(m[j][i]);
        }
      }
      return max(vals);
    }

  template<class utype>
    double frobenius_norm(const xmatrix<utype>& m) {
      double norm = 0;
      for (int i = m.lrows; i <= m.urows; i++) {
        for (int j = m.lcols; j <= m.ucols; j++) {
          norm += abs(m[i][j]) * abs(m[i][j]);
        }
      }
      return sqrt(norm);
    }

  template<class utype>
    double l2_norm(const xmatrix<utype>& m) {
      xvector<utype> wr(m.rows), wi(m.rows);
      eigen(m, wr, wi);
      eigen((trasp(m)*m), wr, wi);
      return sqrt(max(wr));
    }

  template<class utype>
    utype linf_norm(const xmatrix<utype>& m) {
      xvector<utype> vals(m.lrows, m.urows);
      for (int i = m.lrows; i <= m.urows; i++) {
        for (int j = m.lcols; j <= m.ucols; j++) {
          vals[i] += abs(m[i][j]);
        }
      }
      return max(vals);
    }
}

//*****************************************************************************
// ---------------------------------------------------------- aurostd::cematrix
namespace aurostd { // namespace aurostd
#define cematrix_EXIT_RANK_NOT_MATCH 56
#define cematrix_EQUAL_DOUBLE 1.0e-9 // two doubles are equal if difference is smaller than it

  cematrix::cematrix() { // default constructor
    nrow=1;
    ncol=1;
    M=xmatrix<double>(1,1,1,1);
    W=xvector<double>(1,1);
    U=xmatrix<double>(1,1,1,1);
    V=xmatrix<double>(1,1,1,1);
    a_vec=xvector<double>(1,1);
    a_nvec.clear();
    chisq=0.0;
    Cov=xmatrix<double>(1,1,1,1);
  }

  cematrix::cematrix(const xmatrix<double> & A_in) { // copy constructor
    nrow=A_in.rows;
    ncol=A_in.cols;
    M=xmatrix<double>(1,1,nrow,ncol);
    for(int i=1;i<=nrow;i++)
      for(int j=1;j<=ncol;j++)
        M[i][j]=A_in[i][j];
    W=xvector<double>(1,ncol);
    V=xmatrix<double>(1,1,ncol,ncol);
    U=xmatrix<double>(1,1,nrow,ncol);
    a_vec=xvector<double>(1,1);
    a_nvec.clear();
    chisq=0.0;
    Cov=xmatrix<double>(1,1,ncol,ncol);
  }

  cematrix::~cematrix() { // default deconstructor
    a_nvec.clear();
  }

  void cematrix::LeastSquare(xvector<double>& y_vec, xvector<double>& y_sigma) { // function
    if(nrow !=y_vec.rows ) {
      string message = "No match of ranks of b and A. Input two matrices A (m x n) and b (m x 1)";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    //SVDcmp(A);
    SVDFit(y_vec, y_sigma);
  }

  //AS20200811 BEGIN
  void cematrix::LeastSquare(xvector<double> & y_vec){
    xvector<double> y_sigma(y_vec.rows);
    for (int i=y_sigma.lrows; i<=y_sigma.urows; i++) y_sigma[i] = 1.0;
    LeastSquare(y_vec, y_sigma);
  }//AS20200811 END

  double cematrix::Pythag2(double a, double b) { // calculate (a^2+ b^2)^(1/2)
    // from dlapy2.f in Lapack
    double aabs=abs(a),babs=abs(b);
    double val_min,val_max;
    val_min=min(aabs,babs);
    val_max=max(aabs,babs);
    if(val_min ==0) {
      return val_max;
    } else {
      return val_max*sqrt(1.0+(val_min/val_max)*(val_min/val_max));
    }
  }

  void cematrix::SVDsolve(xvector<double>& b_vec) {
    // solve the least squares problem after SVDcmp()
    // it will use xmatrix operators later
    // the results are stored in a_vec
    xvector<double> temp(1,ncol),a_tmp(1,ncol);
    double wj,s;
    const double _NONZERO=1.0e-8;

    for(int j=1;j<=ncol;j++) {
      //wj=W.at(j-1);
      wj=W[j];
      //bj=b_vec[j];
      s=0.0;
      //if(wj !=0.0 )
      if(wj> _NONZERO)
      {
        // only if W[j] !=0
        for(int i=1;i <=nrow;i++) s+=U[i][j]*b_vec[i];
        s /=wj;
        temp[j]=s;
      }
    }

    for(int j=1;j<=ncol;j++) {  //multiply V
      s=0.0;
      for(int i=1;i<=ncol;i++)
        s+=V[j][i]*temp[i];
      a_tmp[j]=s;
    }

    a_vec=a_tmp;
    a_nvec.clear();
    for(int i=1;i<=ncol;i++) {
      a_nvec.push_back(a_vec[i]);
    }
  }

  //void cematrix::SVDcmp_NR()
  bool cematrix::SVDcmp_NR() {
    // SVD decompose matrix A=U Z V^T
    // decomposed matrices U Z V are stored
    xvector<double> rv1(1,ncol);
    double g, scale,anorm;
    int l,nm,jj,j,k,i;
    double f,c,h,x,y,z,s;
    bool flag;
    bool flag_convergence=true;
    // cerr << "ncol " << ncol << " nrow " << nrow << endl;
    // cerr << M << endl;

    l=0;
    g=0.0;
    scale=0.0;
    anorm=0.0;
    xvector<double> W_tmp(1,ncol);
    xmatrix<double> V_tmp(1,1,ncol,ncol);
    xmatrix<double> A(1,1,nrow,ncol);
    for(i=1;i<=nrow;i++)
      for(j=1;j<=ncol;j++)
        if(aurostd::abs(M(i,j))<cematrix_EQUAL_DOUBLE)
          M(i,j)=0;
    A=M; // not destroy input matrix A

    // Householder reduction to bidiagonal form
    for(i=1;i<=ncol;i++) {
      l=i+1;
      rv1[i]=scale*g;
      g=0.0;
      s=0.0;
      scale=0.0;
      if(i <=nrow) {
        for(k=i;k <=nrow;k++ )
          scale+=abs(A[k][i]);
        if(abs(scale)> cematrix_EQUAL_DOUBLE ) { // scale !=0
          for(k=i;k<=nrow;k++) {
            A[k][i] /=scale;
            s+=A[k][i]*A[k][i];
          }
          f=A[i][i];
          g=(-_sign(sqrt(s),f));
          h=f*g - s;
          A[i][i]=f - g;
          for(j=l;j <=ncol;j++) {
            for(s=0.0,k=i;k <=nrow;k++)
              s+=A[k][i]*A[k][j];
            f=s/h;
            for(k=i;k<=nrow;k++)
              A[k][j]+=f*A[k][i];
          }
          for(k=i;k<=nrow ;k++)
            A[k][i] *=scale;
        }
      }
      W_tmp[i]=scale*g;
      g=0.0;
      s=0.0;
      scale=0.0;
      if(i <=nrow && i !=ncol ) {
        for(k=l;k<=ncol;k++)
          scale+=abs(A[i][k]);
        if(scale !=0.0 ) {
          for(k=l;k<=ncol;k++) {
            A[i][k] /=scale;
            s+=A[i][k]*A[i][k];
          }
          f=A[i][l];
          g=(-_sign(sqrt(s),f));
          h=f*g - s;
          A[i][l]=f - g;
          for(k=l;k<=ncol;k++)
            rv1[k]=A[i][k]/h;
          for(j=l;j<=nrow;j++) {
            for(s=0.0,k=l;k <=ncol;k++)
              s+=A[j][k]*A[i][k];
            for(k=l;k<=ncol;k++)
              A[j][k]+=s*rv1[k];
          }
          for(k=l;k<=ncol;k++)
            A[i][k] *=scale;
        }
      }
      anorm=max(anorm,(abs(W_tmp[i])+ abs(rv1[i])) );
    } // i
    for(i=ncol;i>=1;i--) { //Accumulation of right-hand trnasformations
      if(i < ncol) {
        if(g !=0.0 ) {
          for(j=l;j <=ncol;j++ )
            V_tmp[j][i]=(A[i][j]/A[i][l])/g;

          for(j=l;j <=ncol;j++) {
            for(s=0.0,k=l;k<=ncol;k++)
              s+=A[i][k]*V_tmp[k][j];

            for(k=l;k<=ncol;k++)
              V_tmp[k][j]+=s*V_tmp[k][i];

          }
        }
        for(j=l;j<=ncol;j++) {
          V_tmp[i][j]=0.0;
          V_tmp[j][i]=0.0;
        }
      }
      V_tmp[i][i]=1.0;
      g=rv1[i];
      l=i;
    }
    for(i=min(nrow,ncol);i>=1;i--) { // Accumulation of left-hand transformaions
      l=i+1;
      g=W_tmp[i];
      for(j=l;j<=ncol;j++)
        A[i][j]=0.0;

      if(g !=0.0) {
        g=1.0/g;
        for(j=l;j<=ncol;j++) {
          for(s=0.0,k=l;k<=nrow;k++)
            s+=A[k][i]*A[k][j];
          f=(s/A[i][i])*g;
          for(k=i;k<=nrow;k++)
            A[k][j]+=f*A[k][i];
        }
        for(j=i;j<=nrow;j++)
          A[j][i] *=g;
      } else {
        for(j=i;j<=nrow;j++)
          A[j][i]=0.0;
      }
      ++A[i][i];
    }
    for(k=ncol;k>=1;k--) {
      // Diagonalization of the bidiagonal form
      // Loop over singular values and over alowed iterations
      for(int its=1;its<=_MAX_ITS;its++) {
        flag=true;
        //for(l=k;l>=1;l-- ) { // test forsplitting //[CO20200106 - close bracket for indenting]}
        // to avoid out of range as nm cannot be 0
        // should check!!!!!
        nm=1;
        for(l=k;l>=2;l-- ) { // test forsplitting
          nm=l-1;//rv1[1] is always zero
          if((abs(rv1[l])+ anorm)==anorm) {
            flag=false;
            break;
          }
          //if(nm !=0 ) {
          // to avoid abort due to nm ==0
          // xvector start with 1
          if((abs(W_tmp[nm])+ anorm)==anorm ) break;
          //}
        }
        if(flag ) {
          // cancellation of rv1[l],if l> 1
          c=0.0;
          s=1.0;
          for(i=l;i <=k;i++) {
            f=s*rv1[i];
            rv1[i]=c*rv1[i];
            if((abs(f)+ anorm) ==anorm ) break;
            g=W_tmp[i];
            h=Pythag2(f,g);
            W_tmp[i]=h;
            h=1.0/h;
            c=g*h;
            s=-f*h;
            for(j=1;j<=nrow;j++) {
              y=A[j][nm];
              z=A[j][i];
              A[j][nm]=y*c+z*s;
              A[j][i]=z*c-y*s;
            }
          }
        }
        z=W_tmp[k];

        //   cerr << "SVD its=" << its << " l=" << l << " k=" << endl;
        if(l ==k) { // convergence
          //   cerr << "SVD its=" << its << " l=" << l << " k=" << k << endl;
          if(z < 0.0 ) { // singular value is made nonnegative
            W_tmp[k]=-z;
            for(j=1;j<=ncol;j++)
              V_tmp[j][k]=-V_tmp[j][k];
          }
          break;
        }

        //   cerr << "SVD its=" << its << endl;
        if(its ==_MAX_ITS ) {
          cerr << "ERROR - cematrix::SVDcmp_NR: Not converged in " << _MAX_ITS << " SVDcmp iterations !" << endl;;
          cerr << "ERROR - cematrix::SVDcmp_NR: ncol " << ncol << " nrow " << nrow << endl;
          cerr << M << endl;
          flag_convergence=false;
        }
        x=W_tmp[l];// shift from bottom 2x2 minor
        nm=k-1;
        //if(nm !=0 ) {
        y=W_tmp[nm];
        g=rv1[nm];
        //} else {
        //  y=1.0;
        //}
        h=rv1[k];
        f=((y-z)*(y+z)+ (g-h)*(g+h))/(2.0*h*y);
        g=Pythag2(f,1.0);
        f=((x-z)*(x+z)+h*((y/(f+_sign(g,f)))-h))/x;
        c=1.0;s=1.0;// next QR transformation
        for(j=l;j<=nm;j++) {
          i=j+1;g=rv1[i];y=W_tmp[i];h=s*g;g=c*g;
          z=Pythag2(f,h);rv1[j]=z;c=f/z;s=h/z;f=x*c+g*s;
          g=g*c-x*s;h=y*s;y*=c;
          for(jj=1;jj<=ncol;jj++) {
            x=V_tmp[jj][j];
            z=V_tmp[jj][i];
            V_tmp[jj][j]=x*c+z*s;
            V_tmp[jj][i]=z*c-x*s;
          }
          z=Pythag2(f,h);
          W_tmp[j]=z;
          if(z!=0.0) { // rotation can be arbitrary if z=0
            z=1.0/z;c=f*z;s=h*z;
          }
          f=c*g+ s*y;
          x=c*y-s*g;
          for(int jj=1;jj<=nrow;jj++) {
            y=A[jj][j];
            z=A[jj][i];
            A[jj][j]=y*c+z*s;
            A[jj][i]=z*c-y*s;
          }
        }
        rv1[l]=0.0;
        rv1[k]=f;
        W_tmp[k]=x;
      } // its
    } // k

    W=W_tmp;
    V=V_tmp;
    U=A;
    //// output U V W
    //
    //cerr.setf(ios_base::fixed);
    //cerr.precision(6);
    //cerr.width(12);
    ////cout.setw(12);
    //cerr << "*** Decomposition Matrice *** " << endl;;
    //cerr << "*** U matrix *** " << endl;;
    //for(i=1;i <=nrow;i++) { // U
    //  for(j=1;j <=ncol;j++) { cerr << setw(12) << U[i][j] << " "; } cerr << endl;
    //}
    //cerr << "*** W matrix diagonal elements*** " << endl;;
    //for(i=1;i <=ncol;i++) {cerr << setw(12) << W[i] << " ";}
    //cerr << endl;
    //cerr << "*** V matrix *** " << endl;;
    //for(i=1;i <=ncol;i++) { // U
    //  for(j=1;j <=ncol;j++) { cerr << setw(12) << V[i][j] << " ";} cerr << endl;
    //}
    //cerr << "Check the produce against the original matrix " << endl;;
    //cerr << "Original matrix " << endl;;
    //for(i=1;i <=nrow;i++) { // U
    //  for(j=1;j <=ncol;j++) { cerr << setw(12) << M[i][j] << " "; } cerr << endl;
    //}
    //cerr << "Product U W transpose(V) " << endl;;
    for(i=1;i <=nrow;i++) { // U
      for(j=1;j<=ncol;j++) {
        A[i][j]=0.0;
        for(k=1;k <=ncol;k++)
          A[i][j]+=U[i][k]*W[k]*V[j][k];
      }
    }
    //for(i=1;i <=nrow;i++) { // U
    //  for(j=1;j <=ncol;j++) {
    //    cerr << setw(12) << A[i][j] << " ";
    //  }
    //  cerr << endl;
    //}
    //cerr << "A and M is equal " << isequal(A,M) << endl;
    if(isequal(A,M) ==false ) {
      cerr << "ERROR - cematrix::SVDcmp_NR: SVD fails, check the code of cematrix::SVDcmp!" << endl;;
      flag_convergence=false;
    }
    return flag_convergence;

  }

  xmatrix<double> cematrix::InverseMatrix() {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // inverse a general matrix by using SVD
    int i,j,k;

    // SVD to get U,V,W
    SVDcmp_NR();//Two functions are essentially the same
    // inverse of matrix
    xmatrix<double> A_inv(1,1,ncol,nrow);
    double DetW=1.0;// determination of diagonal matrix W
    for(i=1;i <=ncol;i++)
      DetW *=W[i];
    if(DetW < cematrix_EQUAL_DOUBLE) {
      string message = "Singular Matrix. No Inversion.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    for(i=1;i <=ncol;i++) { //row of the inverse Matrix
      for(j=1;j <=nrow;j++) { // colume of the inverse Matrix
        A_inv[i][j]=0.0;
        for(k=1;k <=ncol;k++) // matrix multiplication
          A_inv[i][j]+=V[i][k]/W[k]*U[j][k];
      }
    }
    if(LDEBUG){ //CO20190327
      cerr << "cematrix::InverseMatrix: Inverse of matrix A" << endl;
      for(i=1;i <=ncol;i++) { // U
        for(j=1;j <=nrow;j++)
          cerr << setw(12) << A_inv[i][j] << " ";
        cerr << endl;
      }
    }
    xmatrix<double> Iden_tmp(1,1,nrow,nrow);
    xmatrix<double> Iden(1,1,nrow,nrow);
    for(i=1;i<=nrow;i++) {
      for(j=1;j<=nrow;j++) {
        Iden_tmp[i][j]=0.0;
        if(i ==j ) {
          Iden[i][j]=1.0;
        } else {
          Iden[i][j]=0.0;
        }
        for(k=1;k<=ncol;k++) {
          Iden_tmp[i][j]+=M[i][k]*A_inv[k][j];
        }
      }
    }
    bool verbose = (LDEBUG || isequal(Iden,Iden_tmp));  //CO20190327
    if(verbose){ //CO20190327
      cerr << "cematrix::InverseMatrix: A_inv is the inverse of matrix A? ";
      if(isequal(Iden,Iden_tmp) ) {
        cerr << "Yes!" << endl;;
      } else {
        cerr << "No!" << endl;;
      }
    }
    return A_inv;
  }

  void cematrix::SVDFit(xvector<double>& y,xvector<double>& y_sigma) {
    // Least Square fit by using SVD
    // Here x() is afunc() in numerical recipes in C
    int i,j;
    double wmax,tmp,thresh,sum;
    xmatrix<double> A(1,1,nrow,ncol);
    xmatrix<double> M_orig(1,1,nrow,ncol);
    xvector<double> W_orig(1,ncol);
    xmatrix<double> V_orig(1,1,ncol,ncol);
    xmatrix<double> U_orig(1,1,nrow,ncol);
    //const double TOL=1.0e-13;
    const double TOL=1.0e-8;
    xvector<double> y_cal(1,nrow);
    A=M;
    if(A.rows !=y.rows ) {
      string message = "Ranks of x vector and y vector do not match!";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    xvector<double> b(1,nrow);
    for(i=1;i <=nrow;i++) { // Accumulate coeeficiens of the fitting matrix
      tmp=1.0/y_sigma[i];
      for(j=1;j<=ncol;j++)
        A[i][j]=M[i][j]*tmp;
      b[i]=y[i]*tmp;
    }
    M_orig=M;
    W_orig=W;
    V_orig=V;
    U_orig=U;
    M=A;
    //SVDcmp_NR();// singular value decomposition
    // singular value decomposition
    if(SVDcmp_NR() ) {
      wmax=0.0;
      for(j=1;j<=ncol;j++)
        if(W[j]>wmax) wmax=W[j];
      thresh=TOL*wmax;
      for(j=1;j<ncol;j++)
        if(W[j] < thresh ) W[j]=0.0;
      SVDsolve(b);
      chisq=0.0;
      for(i=1;i<=nrow;i++) {
        for(sum=0.0,j=1;j<=ncol;j++)
          sum+=a_vec[j]*M_orig[i][j];
        tmp=(y[i] - sum)/y_sigma[i];
        chisq+=tmp*tmp;
      }
      chisq=chisq/(nrow - 2.0);// definition in Mathematica
    } else {
      // if svd fails,set the score to a large number
      // to discard it
      chisq=1.0e4;
      for(int i=0;i<y.rows;i++)
        a_nvec.push_back(0.0e0);
    }
    //// output the fitted parameters and square of chi
    //cout << "a_vec " << endl;
    //for(i=1;i<=ncol;i++) {
    //  cout << a_vec[i] << " ";
    //}
    //cout << endl;
    //cout << "chisq " << chisq << endl;
    //cerr << "a_vec " << endl;
    //for(i=1;i<=ncol;i++) {
    //  cerr << a_vec[i] << " ";
    //}
    //cerr << endl;
    //cerr << "chisq " << chisq << endl;
    //cerr << "original y vector" << endl;
    //for(i=1;i<=nrow;i++) {
    //  cerr << y[i] << " ";
    //}
    //  cerr << endl;
    //// fitted value of y vector
    //cerr << "fitted y vector" << endl;
    //for(i=1;i<=nrow;i++) {
    //  for(j=1;j<=ncol;j++) {
    //    y_cal[i]+=A[i][j]*a_vec[j];
    //  }
    //  cerr << y_cal[i] << " ";
    //}
    //cerr << endl;
    //// get the covariance matrix
    //SVDvar();
    //ios_base::fmtflags old_stat=cerr.setf(ios_base::fixed,ios_base::floatfield);
    //cerr.setf(ios_base::scientific);
    //cerr.precision(6);
    //cerr << "Covariance matrix " << endl;;
    //for(i=1;i<=ncol;i++) {
    //  for(j=1;j<=ncol;j++) {
    //    cerr << setw(12) << Cov[i][j] << " ";
    //  }
    //  cerr << endl;
    //}

    // set everything back to original values
    M=M_orig;
    U=U_orig;
    V=V_orig;
    W=W_orig;
    //cerr.setf(old_stat,ios_base::floatfield);
    //cerr.unsetf(ios_base::scientific);

  }

  void cematrix::SVDvar() {
    // get the covariance matrix Cov
    xmatrix<double> Cov_tmp(1,1,ncol,ncol);
    int k,j,i;
    double sum;
    xvector<double> wti(1,ncol);
    for(i=1;i<= ncol;i++) {
      wti[i]=0.0;
      if(W[i] !=0.0 )
        wti[i]=1.0/(W[i]*W[i]);
    }
    for(i=1;i<=ncol;i++) { // sum contributions to covariance matrix
      for(j=1;j<=i;j++) {
        for(sum=0.0,k=1;k<=ncol;k++)
          sum+=V[i][k]*V[j][k]*wti[k];
        Cov_tmp[j][i]=sum;
        Cov_tmp[i][j]=sum;
      }
    }
    Cov=Cov_tmp*chisq;// definition in Mathematica
  }

  xvector<double> cematrix::EigenValues() {
    // inverse a general matrix by using SVD
    // SVD to get U,V,W
    SVDcmp_NR(); //Two functions are essentially the same
    return W;
  }
}

// **************************************************************************
// **************************************************************************
//CO20200404 - moving matrix() from pflow to aurostd because it is templated
//doesn't compile otherwise

// ***************************************************************************
// Matrix classes in Dane Morgan Style
// ***************************************************************************
namespace aurostd {
  // constructor

  template<class utype> matrix<utype>::matrix(void) {free();} // default

  template<class utype>
    matrix<utype>::matrix(const int m) { // specifying rows
      free();  //CO20200404 pflow::matrix()->aurostd::matrix()
      std::vector<utype> v;
      mat=std::vector<std::vector<utype> > (m,v);
    }

  template<class utype>
    matrix<utype>::matrix(const int m, const int n) { // specifying rows and columns
      free();  //CO20200404 pflow::matrix()->aurostd::matrix()
      std::vector<utype> v(n);
      mat=std::vector<std::vector<utype> > (m,v);
    }

  template<class utype>
    matrix<utype>::matrix(const int m, const std::vector<utype>& inutypevec) { // specifying rows as vectors.
      free();  //CO20200404 pflow::matrix()->aurostd::matrix()
      mat=std::vector<std::vector<utype> > (m,inutypevec);
    }

  template<class utype>
    matrix<utype>::matrix(const int m, const int n, const utype& inutype) { // specifying rows and columns and initial values
      free();  //CO20200404 pflow::matrix()->aurostd::matrix()
      std::vector<utype> v(n);
      mat=std::vector<std::vector<utype> > (m,v);
      for(uint i=0;i<mat.size();i++) {
        for(uint j=0;j<mat[i].size();j++) {
          mat[i][j]=inutype;
        }
      }
    }

  //CO20200404 START - patching matrix for nietzsche
  template<class utype> matrix<utype>::matrix(const matrix& b){copy(b);}
  template<class utype> matrix<utype>::~matrix(void){free();}
  template<class utype> void matrix<utype>::clear() {matrix a;copy(a);}  //clear PUBLIC

  template<class utype>
    const matrix<utype>& matrix<utype>::operator=(const matrix& other) {
      if(this!=&other) {copy(other);}
      return *this;
    }

  template<class utype>
    void matrix<utype>::free(){
      for(uint i=0;i<mat.size();i++){mat[i].clear();} mat.clear();
    }

  template<class utype>
    void matrix<utype>::copy(const matrix& b){
      free();
      for(uint i=0;i<b.mat.size();i++){mat.push_back(std::vector<utype>(0));for(uint j=0;j<b.mat[i].size();j++){mat[i].push_back(b.mat[i][j]);}}
    }
  //CO20200404 STOP - patching matrix for nietzsche

  // accessors
  template<class utype>
    void matrix<utype>::print(void) {
      cout.setf(std::ios::fixed,std::ios::floatfield);
      cout.precision(4);
      for(uint i=0;i<mat.size();i++) {
        cout << "  ";
        for(uint j=0;j<mat[i].size();j++) {
          cout << " " << mat[i][j];
        }
        cout << endl;
      }
    }

  //   template<class utype>
  //   inline int matrix<utype>::size(void) const{
  //     return (int) mat.size();
  //   }


  template<class utype>
    matrix<utype> matrix<utype>::transpose() const{
      matrix<utype> tmat;
      uint m=mat.size();
      if(m==0) return tmat;
      uint n=mat[0].size();
      tmat = matrix<utype> (n,m);
      for(uint i=0;i<m;i++) {
        for(uint j=0;j<n;j++) {
          tmat[j][i]=mat[i][j];
        }
      }
      return tmat;
    }

  // template<class utype>
  // std::vector<std::vector<utype> >::iterator matrix<utype>::begin() {
  // return mat.begin();
  // }
  //
  // template<class utype>
  // std::vector<std::vector<utype> >::iterator matrix<utype>::end() {
  // return mat.end();
  // }

  // operator
  //   template<class utype>
  //   std::vector<utype>& matrix<utype>::operator[] (const int index) {
  //     assert(index>=0 && index<=mat.size());
  //     return mat[index];
  //   }
  //   template<class utype>
  //   const std::vector<utype>& matrix<utype>::operator[] (const int index) const {
  //     assert(index>=0 && index<=mat.size());
  //     return mat[index];

  //   }

  //[CO20200404 - OBSOLETE]template<class utype>
  //[CO20200404 - OBSOLETE]  const matrix<utype>& matrix<utype>::operator=(const matrix<utype> &b) {
  //[CO20200404 - OBSOLETE]    if(this != &b) {
  //[CO20200404 - OBSOLETE]      uint m=b.mat.size();
  //[CO20200404 - OBSOLETE]      uint n=0;
  //[CO20200404 - OBSOLETE]      mat=std::vector<std::vector<utype> > (m);
  //[CO20200404 - OBSOLETE]      for(uint i=0;i<m;i++) {
  //[CO20200404 - OBSOLETE]        n=b.mat[i].size();
  //[CO20200404 - OBSOLETE]        mat[i]=std::vector<utype> (n);
  //[CO20200404 - OBSOLETE]        for(uint j=0;j<n;j++) {
  //[CO20200404 - OBSOLETE]          mat[i][j]=b.mat[i][j];
  //[CO20200404 - OBSOLETE]        }
  //[CO20200404 - OBSOLETE]      }
  //[CO20200404 - OBSOLETE]    }
  //[CO20200404 - OBSOLETE]    return *this;
  //[CO20200404 - OBSOLETE]  }

  // mutators
  //  template<class utype>
  // void matrix<utype>::push_back(const std::vector<utype>& inutypevec) {
  //   mat.push_back(inutypevec);
  //  }

  //  template<class utype>
  // void matrix<utype>::pop_back() {
  //   mat.pop_back();
  //  }

  template<class utype>
    void matrix<utype>::vecvec2mat(const std::vector<std::vector<utype> >& inVV) {
      mat=std::vector<std::vector<utype> > (inVV.size());
      for(uint i=0;i<mat.size();i++) {
        mat[i]=std::vector<utype> (inVV[i].size());
        for(uint j=0;j<mat[i].size();j++) {
          mat[i][j]=inVV[i][j];
        }
      }
    }

  // template<class utype>
  // void matrix<utype>::clear() {
  //   mat.clear();
  // }

  template<class utype>
    void matrix<utype>::vec2mat(const std::vector<utype>& inV) {
      mat=std::vector<std::vector<utype> > (1);
      mat[0]=std::vector<utype> (inV.size());
      for(uint j=0;j<mat[0].size();j++) {
        mat[0][j]=inV[j];
      }
    }

  // template<class utype>
  // void matrix<utype>::insert(const int& id, const std::vector<utype>& inV) {
  //  mat.insert(mat.begin()+id,inV);
  // }

  // template<class utype>
  // void matrix<utype>::erase(const int id) {
  // std::vector<utype>::iterator p=mat.begin()+id;
  // mat.erase(p);
  // }

  // template<class utype>
  // void matrix<utype>::erase_col(const int id) {
  //   for(int i=0;i<mat.size();i++) {
  //           std::vector<utype>::iterator p=mat[i].begin()+id;
  //     if(id<mat[i].size()) mat[i].erase(p);
  //   }
  // }

  // template<class utype>
  // void matrix<utype>::erase(const int id) {
  //   std::vector<std::vector<utype> >::iterator p=mat.begin()+id;
  //   mat.erase(p);
  // }

  template <class utype> matrix<utype>
    xmatrix2matrix(const xmatrix<utype>& _xmatrix) {
      int isize=_xmatrix.rows,jsize=_xmatrix.cols;
      matrix<utype> _matrix(isize,jsize);
      for(int i=0;i<isize;i++)   //HE20220124 removed register as it is deprecated in C++11 and gone in C++17
        for(int j=0;j<jsize;j++) //Defect report 809 http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n4193.html#809
          _matrix[i][j]=_xmatrix(i+_xmatrix.lrows,j+_xmatrix.lcols);
      return _matrix;
    }

  //   matrix<double> xmatrix2matrix(const xmatrix<double>& _xmatrix) {
  //     int isize=_xmatrix.rows,jsize=_xmatrix.cols;
  //     matrix<double> _matrix(isize,jsize);
  //     for(int i=0;i<isize;i++)
  //       for(int j=0;j<jsize;j++)
  // 	_matrix[i][j]=_xmatrix(i+_xmatrix.lrows,j+_xmatrix.lcols);
  //     return _matrix;
  //   }

  template <class utype> xmatrix<utype>
    matrix2xmatrix(const matrix<utype>& _matrix) {
      int isize=_matrix.size(),jsize=_matrix[0].size();
      xmatrix<utype> _xmatrix(isize,jsize);
      for(int i=1;i<=isize;i++)   //HE20220124 removed register as it is deprecated in C++11 and gone in C++17
        for(int j=1;j<=jsize;j++) //Defect report 809 http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n4193.html#809
          _xmatrix(i,j)=_matrix[i-1][j-1];
      return _xmatrix;
    }



}

// **************************************************************************
namespace aurostd { //force the compiler to instantiate the template at this point (avoids linker issues)
  template class xmatrix<int>;
  template class xmatrix<unsigned int>;
  template class xmatrix<long int>;
  template class xmatrix<long unsigned int>;
  template class xmatrix<long long int>;
  template class xmatrix<long long unsigned int>;
  template class xmatrix<float>;
  template class xmatrix<double>;
  template class xmatrix<long double>;
}


#endif
// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
