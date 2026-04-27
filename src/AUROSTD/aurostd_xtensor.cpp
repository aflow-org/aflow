//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************
//
// This class implements n-dimensional tensors using a 1D array to mimic a
// tensor. See aurostd_xtensor.h for descriptions of the class attributes.

#include "aurostd_xtensor.h"

#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "aurostd.h"
#include "aurostd_automatic_template.h"
#include "aurostd_xerror.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xscalar.h"
#include "aurostd_xvector.h"

#include "aflow_xhost.h" // todo required for XHOST.DEBUG use

#define _DEBUG_XTENSOR_ false

// _CHECK_BOUNDS_ defines whether a bounds check is performed when indexing an
// xtensor. While it is safe to set it to true, it greatly impact performance
// It is recommended to uncomment the line below when testing new code and set
// comment it out once the code it working.
// #define _CHECK_BOUNDS_

static const std::string _SUBTENSOR_ERR_PREFIX_ = "aurostd::_subtensor::";
static const std::string _XTENSOR_ERR_PREFIX_ = "aurostd::xtensor::";

using std::function;
using std::ifstream;
using std::iostream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;

using aurostd::xerror;

/********************************* SUBTENSOR ********************************/
// Helper class for indexing of xtensor. Indexing xtensor is performed by
// storing the indices subsequently into _subtensor and then by using the
// assignment operator or type casting (see conversion operators).
//
// The indices in _subtensor must be 0-based, so they need to be converted
// before being passed into the constructor or appended to the indices vector.
namespace aurostd {

  // Constructors////////////////////////////////////////////////////////////////
  template <class utype> _subtensor<utype>::_subtensor(const int& i, utype* crp, const xtensor<utype>& tnsr) : _tensor(tnsr) {
    indexed_dim = 1;
    shift = (i - _tensor.lindex[0]) * _tensor.shifts[0];
    corpus = crp;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>::_subtensor(const int&, utype*, const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>::_subtensor(const std::vector<int>& ind, utype* crp, const xtensor<utype>& tnsr) : _tensor(tnsr) {
    const uint isize = ind.size();
    indexed_dim = isize;
    shift = 0;
    for (uint i = 0; i < isize; i++) {
      shift += (ind[i] - _tensor.lindex[i]) * _tensor.shifts[i];
    }
    corpus = crp;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>::_subtensor(const std::vector<int>&, utype*, const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>::_subtensor(const aurostd::xvector<int>& ind, utype* crp, const xtensor<utype>& tnsr) : _tensor(tnsr) {
    indexed_dim = ind.rows;
    shift = 0;
    for (int i = 0, ixvec = ind.lrows; ixvec <= ind.urows; i++, ixvec++) {
      shift += (ind[ixvec] - _tensor.lindex[i]) * _tensor.shifts[i];
    }
    corpus = crp;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>::_subtensor(const aurostd::xvector<int>&, utype*, const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Destructor//////////////////////////////////////////////////////////////////
  template <class utype> _subtensor<utype>::~_subtensor() {}
#define AST_TEMPLATE(utype) template _subtensor<utype>::~_subtensor();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Indexing operator///////////////////////////////////////////////////////////
  template <class utype> _subtensor<utype>& _subtensor<utype>::operator[](const int& i) {
#ifdef _CHECK_BOUNDS_
    std::stringstream message;
    if (indexed_dim == _tensor.ndim) {
      message << "Cannot subscribe tensor any further (tensor size: " << _tensor.ndim << ").";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    } else if ((i < _tensor.lindex[indexed_dim]) || (i > _tensor.uindex[indexed_dim])) {
      message << "Index " << i << " out of bounds for dimension " << (indexed_dim);
      message << " (lindex = " << _tensor.lindex[indexed_dim] << ", uindex = ";
      message << _tensor.uindex[indexed_dim] << ")";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    shift += (i - _tensor.lindex[indexed_dim]) * _tensor.shifts[indexed_dim];
    indexed_dim++;
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator[](const int&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator()(const std::vector<int>& ind) {
    const uint isize = ind.size();
#ifdef _CHECK_BOUNDS_
    std::stringstream message;
    if (indexed_dim + isize > _tensor.ndim) {
      message << "Too many indices for tensor with " << _tensor.ndim << " dimensions";
      message << " (" << (indexed_dim + isize) << " provided).";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
#endif
    for (uint i = 0; i < isize; i++) {
#ifdef _CHECK_BOUNDS_
      if ((ind[i] < _tensor.lindex[i + indexed_dim]) || (ind[i] > _tensor.uindex[i + indexed_dim])) {
        message << "Index " << ind[i] << " out of bounds for dimension " << (i + indexed_dim);
        message << " (lindex = " << _tensor.lindex[i + indexed_dim];
        message << ", uindex = " << _tensor.uindex[i + indexed_dim] << ")";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
#endif
      shift += (ind[i] - _tensor.lindex[indexed_dim]) * _tensor.shifts[indexed_dim];
      indexed_dim++;
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator()(const std::vector<int>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator()(const aurostd::xvector<int>& ind) {
#ifdef _CHECK_BOUNDS_
    std::stringstream message;
    if (indexed_dim + ind.rows > _tensor.ndim) {
      message << "Too many indices for tensor with " << _tensor.ndim << " dimensions";
      message << " (" << (indexed_dim + ind.rows) << " provided).";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
#endif
    for (int i = ind.lrows; i <= ind.urows; i++) {
#ifdef _CHECK_BOUNDS_
      if ((ind[i] < _tensor.lindex[indexed_dim]) || (ind[i] > _tensor.uindex[indexed_dim])) {
        message << "Index " << ind[i] << " out of bounds for dimension " << indexed_dim;
        message << " (lindex = " << _tensor.lindex[indexed_dim];
        message << ", uindex = " << _tensor.uindex[indexed_dim] << ")";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
#endif
      shift += (ind[i] - _tensor.lindex[indexed_dim]) * _tensor.shifts[indexed_dim];
      indexed_dim++;
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator()(const aurostd::xvector<int>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void _subtensor<utype>::set(const utype& val) {
    const int n = getNumElements();
    for (int s = 0; s < n; s++) {
      corpus[shift + s] = val;
    }
  }
#define AST_TEMPLATE(utype) template void _subtensor<utype>::set(const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void _subtensor<utype>::set(const utype& val, const int& i) {
    corpus[i] = val;
  }
#define AST_TEMPLATE(utype) template void _subtensor<utype>::set(const utype&, const int&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype _subtensor<utype>::get(const int& i) const {
    return corpus[i];
  }
#define AST_TEMPLATE(utype) template utype _subtensor<utype>::get(const int&) const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // BEGIN Operators
  // Unary operators/////////////////////////////////////////////////////////////
  //// BEGIN with scalars
  template <class utype> _subtensor<utype>& _subtensor<utype>::operator+=(utype scalar) {
    if (indexed_dim == _tensor.ndim) {
      corpus[shift] += scalar;
    } else {
      const int n = getNumElements();
      for (int s = 0; s < n; s++) {
        corpus[shift + s] += scalar;
      }
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator+=(utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator-=(utype scalar) {
    if (indexed_dim == _tensor.ndim) {
      corpus[shift] -= scalar;
    } else {
      const int n = getNumElements();
      for (int s = 0; s < n; s++) {
        corpus[shift + s] -= scalar;
      }
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator-=(utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator*=(utype scalar) {
    if (indexed_dim == _tensor.ndim) {
      corpus[shift] *= scalar;
    } else {
      const int n = getNumElements();
      for (int s = 0; s < n; s++) {
        corpus[shift + s] *= scalar;
      }
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator*=(utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator/=(utype scalar) {
    if (indexed_dim == _tensor.ndim) {
      corpus[shift] /= scalar;
    } else {
      const int n = getNumElements();
      for (int s = 0; s < n; s++) {
        corpus[shift + s] /= scalar;
      }
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator/=(utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  //// END with scalars

  //// BEGIN with tensors

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator+=(const _subtensor<utype>& st) {
    if (sameShape(st)) {
      if (indexed_dim == _tensor.ndim) {
        corpus[shift] += st.get(st.shift);
      } else {
        const int n = getNumElements();
        for (int s = 0; s < n; s++) {
          corpus[shift + s] += st.get(st.shift + s);
        }
      }
    } else {
      const std::string message = "Subtensors have different shapes.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator+=(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator-=(const _subtensor<utype>& st) {
    if (sameShape(st)) {
      if (indexed_dim == _tensor.ndim) {
        corpus[shift] -= st.corpus[st.shift];
      } else {
        const int n = getNumElements();
        for (int s = 0; s < n; s++) {
          corpus[shift + s] -= st.get(st.shift + s);
        }
      }
    } else {
      const std::string message = "Subtensors have different shapes.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator-=(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator+=(const xtensor<utype>& tensor) {
    if (sameShape(tensor)) {
      if (indexed_dim == _tensor.ndim) {
        corpus[shift] += tensor[0];
      } else {
        for (int i = 0; i < tensor.nelements; i++) {
          corpus[shift + i] += tensor.get(i);
        }
      }
    } else {
      const std::string message = "Subtensor and tensor have different shapes.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator+=(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype>& _subtensor<utype>::operator-=(const xtensor<utype>& tensor) {
    if (sameShape(tensor)) {
      if (indexed_dim == _tensor.ndim) {
        corpus[shift] -= tensor[0];
      } else {
        for (int i = 0; i < tensor.nelements; i++) {
          corpus[shift + i] -= tensor.get(i);
        }
      }
    } else {
      const std::string message = "Subtensor and tensor have different shapes.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator-=(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  //// END with tensors

  //// BEGIN with subtensors

  template <class utype> xtensor<utype> operator+(const _subtensor<utype>& st) {
    return xtensor<utype>(st);
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator+(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator-(const _subtensor<utype>& st) {
    const xtensor<utype> tensor(st);
    return -tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator-(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  //// END with subtensors

  // Conversion operators////////////////////////////////////////////////////////
  //  Converts a _subtensor into utype to be assigned to a utype variable. For
  //  the conversion of _subtensor into xtensor, see the (copy) constructor for
  //  xtensor.
  template <class utype> _subtensor<utype>::operator utype() const {
    return corpus[shift];
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>::operator utype() const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Assigns a utype to an xtensor element
  template <class utype> _subtensor<utype>& _subtensor<utype>::operator=(const utype& val) {
    if (indexed_dim == _tensor.ndim) {
      set(val, shift);
        // Exception handling
    } else {
      std::stringstream message;
      message << "Cannot assign single value to tensor of size " << (_tensor.ndim - indexed_dim) << ".";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator=(const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Assigns a subtensor to another subtensor
  template <class utype> _subtensor<utype>& _subtensor<utype>::operator=(const _subtensor<utype>& st) {
    if (sameShape(st)) {
      const int end = _tensor.shifts[indexed_dim - 1] + shift;
      for (int i = shift, ist = st.shift; i < end; i++, ist++) {
        corpus[i] = st.get(ist);
      }
    } else {
      const std::string message = "Subtensors have different shapes.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator=(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Assigns a tensor to a subtensor
  template <class utype> _subtensor<utype>& _subtensor<utype>::operator=(const xtensor<utype>& tensor) {
    if (sameShape(tensor)) {
      for (int i = 0, ist = shift; i < tensor.nelements; i++, ist++) {
        corpus[ist] = tensor.get(i);
      }
    } else {
      const std::string message = "Subtensor and tensor have different shapes.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator=(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Assigns a vector to a 1D slice of a tensor
  template <class utype> _subtensor<utype>& _subtensor<utype>::operator=(const std::vector<utype>& vec) {
    bool throw_exception = false;
    std::stringstream message;
      // Slice must be 1D
    if (indexed_dim == (_tensor.ndim - 1)) {
        // Tensor slice and vector must have same size
      const int vsize = (int) vec.size();
      if (vsize == _tensor.shape[indexed_dim]) {
        for (int i = 0, ist = shift; i < vsize; i++, ist++) {
          corpus[ist] = vec[i];
        }
          // Exception handling.
      } else {
        throw_exception = true;
        message << "Tensor slice and vector have different sizes. ";
        message << "Size of tensor slice: " << _tensor.shape[indexed_dim] << "; ";
        message << "size of vector: " << vec.size() << ".";
      }
    } else {
      throw_exception = true;
      message << "Tensor slice must be 1D to assign vector (is " << (_tensor.ndim - indexed_dim) << "D).";
    }
    if (throw_exception) {
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator=(const std::vector<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Assigns an xvector to a 1D slice of a tensor
  template <class utype> _subtensor<utype>& _subtensor<utype>::operator=(const aurostd::xvector<utype>& xvec) {
    bool throw_exception = false;
    std::stringstream message;
      // Tensor slice must be 1D
    if (indexed_dim == (_tensor.ndim - 1)) {
        // Tensor slice and xvector must have same size
      if (xvec.rows == _tensor.shape[indexed_dim]) {
        for (int i = xvec.lrows, ist = shift; i <= xvec.urows; i++, ist++) {
          corpus[ist] = xvec[i];
        }
          // Exception handling
      } else {
        throw_exception = true;
        message << "Tensor slice and xvector have different sizes. ";
        message << "Size of tensor slice: " << _tensor.shape[indexed_dim] << "; ";
        message << "size of xvector: " << xvec.rows << ".";
      }
    } else {
      throw_exception = true;
      message << "Tensor slice must be 1D to assign xvector (is " << (_tensor.ndim - indexed_dim) << "D).";
    }
    if (throw_exception) {
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator=(const aurostd::xvector<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Assigns an xmatrix to a 2D slice of a tensor
  template <class utype> _subtensor<utype>& _subtensor<utype>::operator=(const aurostd::xmatrix<utype>& xmat) {
    bool throw_exception = false;
    std::stringstream message;
      // Tensor slice must be 2D
    if (indexed_dim == (_tensor.ndim - 2)) {
        // Tensor slice must have same shape as matrix
      if ((xmat.rows == _tensor.shape[indexed_dim]) && (xmat.cols == _tensor.shape[indexed_dim + 1])) {
        for (int i = xmat.lrows; i <= xmat.urows; i++) {
          for (int j = xmat.lcols; j <= xmat.ucols; j++) {
            corpus[shift] = xmat[i][j];
            shift++;
          }
        }
      } else {
        throw_exception = true;
        message << "Tensor slice and xmatrix have different shapes. ";
        message << "Shape of tensor slice: (" << _tensor.shape[indexed_dim] << ", " << _tensor.shape[indexed_dim + 1] << "); ";
        message << "shape of xmatrix: (" << xmat.rows << ", " << xmat.cols << ").";
      }
    } else {
      throw_exception = true;
      message << "Tensor slice must be 2D to assign xmatrix (is " << (_tensor.ndim - indexed_dim) << "D).";
    }
    if (throw_exception) {
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template _subtensor<utype>& _subtensor<utype>::operator=(const aurostd::xmatrix<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  // END Operators

  // BEGIN Helper functions
  // getNumElements//////////////////////////////////////////////////////////////
  template <class utype> int _subtensor<utype>::getNumElements() const {
    int n = 1;
    for (uint i = indexed_dim; i < _tensor.ndim; i++) {
      n *= _tensor.shape[i];
    }
    return n;
  }
#define AST_TEMPLATE(utype) template int _subtensor<utype>::getNumElements() const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // sameShape///////////////////////////////////////////////////////////////////
  template <class utype> bool _subtensor<utype>::sameShape(const _subtensor<utype>& st) {
    const uint size1 = _tensor.ndim - indexed_dim;
    const uint size2 = st._tensor.ndim - st.indexed_dim;
    if (size1 != size2) {
      return false;
    } else {
      const int diff = (int) (size2 - size1);
      for (uint i = indexed_dim; i < _tensor.ndim; i++) {
        if (_tensor.shape[i] != st._tensor.shape[i + diff]) {
          return false;
        }
      }
    }
    return true;
  }
#define AST_TEMPLATE(utype) template bool _subtensor<utype>::sameShape(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> bool _subtensor<utype>::sameShape(const xtensor<utype>& tensor) {
    const uint size = _tensor.ndim - indexed_dim;
    if (size != tensor.ndim) {
      return false;
    } else {
      for (uint i = 0; i < size; i++) {
        if (_tensor.shape[i + indexed_dim] != tensor.shape[i]) {
          return false;
        }
      }
    }
    return true;
  }
#define AST_TEMPLATE(utype) template bool _subtensor<utype>::sameShape(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // END Helper functions
} // namespace aurostd

/******************************* END SUBTENSOR ******************************/

/********************************** XTENSOR *********************************/
// The xtensor class. As with xvector and xmatrix, xtensors are initialized
// upper indices (uind/uindex) and the lower indices (lind/lindex), and the
// default indexing is 1-based. uindex and lindex can be provided using
// (x)vectors or an integer describing the number of dimensions and the lower
// index for all dimensions (default 1).
namespace aurostd {

  // BEGIN Constructors/Destructors
  // Constructors////////////////////////////////////////////////////////////////
  template <class utype> xtensor<utype>::xtensor() {  // Default constructor
    const std::vector<int> uind(3, 3);
    const std::vector<int> lind(3, 1);
    buildTensor(uind, lind);
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::xtensor();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>::xtensor(const std::vector<int>& uind) {
    const std::vector<int> lind(uind.size(), 1);
    buildTensor(uind, lind);
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::xtensor(const std::vector<int>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>::xtensor(const std::vector<int>& uind, const std::vector<int>& lind) {
    buildTensor(uind, lind);
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::xtensor(const std::vector<int>&, const std::vector<int>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>::xtensor(const aurostd::xvector<int>& xuind) {
    const std::vector<int> uind = aurostd::xvector2vector(xuind);
    const std::vector<int> lind(uind.size(), 1);
    buildTensor(uind, lind);
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::xtensor(const aurostd::xvector<int>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>::xtensor(const aurostd::xvector<int>& xuind, const aurostd::xvector<int>& xlind) {
    const std::vector<int> uind = aurostd::xvector2vector(xuind);
    const std::vector<int> lind = aurostd::xvector2vector(xlind);
    buildTensor(uind, lind);
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::xtensor(const aurostd::xvector<int>&, const aurostd::xvector<int>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>::xtensor(const int& ui, const int& li) {
      // Creates a tensor with (ui-li + 1) dimensions of size (ui-li). For example,
      // a 3x3x3 tensor would be xtensor(3) (1-based) or xtensor(2, 0) (0-based).
    const int n = ui - li + 1;
    const std::vector<int> uind(n, n);
    const std::vector<int> lind(n, li);
    buildTensor(uind, lind);
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::xtensor(const int&, const int&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>::xtensor(const _subtensor<utype>& st) {
    std::vector<int> uind;
    std::vector<int> lind;
    if (st.indexed_dim == st._tensor.ndim) {
      uind.resize(1, 1);
      lind.resize(1, st._tensor.lindex[st._tensor.ndim - 1]);
    } else {
      const int size = (int) st._tensor.ndim - st.indexed_dim;
      uind.resize(size);
      lind.resize(size);
      for (int i = 0; i < size; i++) {
        uind[i] = st._tensor.uindex[i + st.indexed_dim];
        lind[i] = st._tensor.lindex[i + st.indexed_dim];
      }
    }
    buildTensor(uind, lind);
    const int end = st.shift + st._tensor.shifts[st.indexed_dim - 1];
    for (int i = 0, ist = st.shift; ist < end; i++, ist++) {
      corpus[i] = st.get(ist);
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::xtensor(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // buildTensor/////////////////////////////////////////////////////////////////
  template <class utype> void xtensor<utype>::buildTensor(const std::vector<int>& uind, const std::vector<int>& lind) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (checkInit(uind, lind)) {
      ndim = uind.size();
      uindex = new int[ndim];
      lindex = new int[ndim];
      shape = new int[ndim];
      shifts = new int[ndim];

      for (uint i = 0; i < ndim; i++) {
        uindex[i] = uind[i];
        lindex[i] = lind[i];
        shape[i] = uindex[i] - lindex[i] + 1;
        shifts[i] = 1;
      }

      if (LDEBUG) {
        std::cerr << _XTENSOR_ERR_PREFIX_ << "buildTensor - Tensor dimensions:";
        for (uint i = 0; i < ndim; i++) {
          std::cerr << " " << shape[i];
        }
        std::cerr << std::endl;
      }

      for (uint s = 0; s < ndim - 1; s++) {
        for (uint i = s + 1; i < ndim; i++) {
          shifts[s] *= shape[i];
        }
      }

        // Special case for a 0D tensor
      if ((uind.size() == 1) && ((uindex[0] - lindex[0]) == 0)) {
        ndim = 0;
      }

      if (ndim < 2) {
        is_cubic = false;
      } else {
        is_cubic = true;
        for (uint i = 1; i < ndim; i++) {
          if (shape[i] != shape[i - 1]) {
            is_cubic = false;
            i = ndim;
          }
        }
      }

      nelements = 1;
      for (uint i = 0; i < ndim; i++) {
        nelements *= shape[i];
      }

      size = (char) (sizeof(utype));
      tsize = (long int) size * nelements;

      corpus = new utype[nelements];
      for (int i = 0; i < nelements; i++) {
        corpus[i] = (utype) 0;
      }

      if (LDEBUG) {
        std::cerr << "is_cubic = " << ((is_cubic) ? "true" : "false");
        std::cerr << ", sizeof = " << size << ", tsize = " << tsize << std::endl;
      }
      if (!corpus) {
        const std::string message = "Allocation failure in xtensor constructor.";
        throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _ALLOC_ALLOCATE_);
      }
        // Exception handling
    } else {
      tsize = 0;  // For destructor
      std::stringstream message;
      int code;
      if (uind.size() != lind.size()) {
        message << "uindex and lindex are of different size.";
        code = _INDEX_MISMATCH_;
      } else {
        if (uind.empty() || lind.empty()) {
          message << "Cannot initialize tensor with empty vector.";
        } else {
          message << "Vector describing tensor indices contains illegal values (uindex < lindex).";
        }
        code = _INDEX_ILLEGAL_;
      }
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, code);
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::buildTensor(const std::vector<int>&, const std::vector<int>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Copy Constructors///////////////////////////////////////////////////////////
  template <class utype> xtensor<utype>::xtensor(const xtensor<utype>& that) {
    is_cubic = that.is_cubic;
    nelements = that.nelements;
    ndim = that.ndim;
    tsize = that.tsize;

    uindex = new int[ndim];
    lindex = new int[ndim];
    shape = new int[ndim];
    shifts = new int[ndim];
    for (uint i = 0; i < ndim; i++) {
      uindex[i] = that.uindex[i];
      lindex[i] = that.lindex[i];
      shape[i] = that.shape[i];
      shifts[i] = that.shifts[i];
    }

    corpus = new utype[nelements];
    for (int i = 0; i < nelements; i++) {
      corpus[i] = that.corpus[i];
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::xtensor(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> xtensor<utype>::operator=(const xtensor<utype>& that) {
    if (corpus != that.corpus) {
      is_cubic = that.is_cubic;
      ndim = that.ndim;
      nelements = that.nelements;
      size = that.size;
      tsize = that.tsize;

      delete[] (uindex);
      delete[] (lindex);
      delete[] (shape);
      delete[] (shifts);
      uindex = new int[ndim];
      lindex = new int[ndim];
      shape = new int[ndim];
      shifts = new int[ndim];
      for (uint i = 0; i < ndim; i++) {
        uindex[i] = that.uindex[i];
        lindex[i] = that.lindex[i];
        shape[i] = that.shape[i];
        shifts[i] = that.shifts[i];
      }

      delete[] (corpus);
      corpus = new utype[nelements];
      for (int i = 0; i < nelements; i++) {
        corpus[i] = that.corpus[i];
      }
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> xtensor<utype>::operator=(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> xtensor<utype>::operator=(const _subtensor<utype>& st) {
    const xtensor<utype> tensor(st);
    *this = tensor;
    return *this;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> xtensor<utype>::operator=(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // Destructor//////////////////////////////////////////////////////////////////
  template <class utype> xtensor<utype>::~xtensor() {
    if (tsize > 0) {
      delete[] (corpus);
      delete[] (shape);
      delete[] (uindex);
      delete[] (lindex);
      delete[] (shifts);
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype>::~xtensor();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  // END Constructors/Destructors

  // BEGIN helper functions
  // checkInit///////////////////////////////////////////////////////////////////
  // helper function for the constructor to test the validity of the input std::vector
  template <class utype> bool xtensor<utype>::checkInit(const std::vector<int>& uind, const std::vector<int>& lind) {
    if (uind.size() != lind.size()) {
      return false;
    }
    if (uind.empty() || lind.empty()) {
      return false;
    }
    for (size_t i = 0; i < uind.size(); i++) {
      if (uind[i] < lind[i]) {
        return false;
      }
    }
    return true;
  }
#define AST_TEMPLATE(utype) template bool xtensor<utype>::checkInit(const std::vector<int>&, const std::vector<int>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // sameShape///////////////////////////////////////////////////////////////////
  //  helper function for operations between two tensors to check if both have
  //  the same shape
  template <class utype> bool xtensor<utype>::sameShape(const xtensor<utype>& tensor) const {
    for (uint d = 0; d < ndim; d++) {
      if (shape[d] != tensor.shape[d]) {
        return false;
      }
    }
    return true;
  }
#define AST_TEMPLATE(utype) template bool xtensor<utype>::sameShape(const xtensor<utype>&) const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> bool xtensor<utype>::sameShape(const _subtensor<utype>& st) const {
    const uint size = st._tensor.ndim - (uint) st.indexed_dim;
    if (size != ndim) {
      return false;
    } else {
      for (uint i = 0; i < size; i++) {
        if (st._tensor.shape[i + st.indexed_dim] != shape[i]) {
          return false;
        }
      }
    }
    return true;
  }
#define AST_TEMPLATE(utype) template bool xtensor<utype>::sameShape(const _subtensor<utype>&) const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  // END helper functions

  // BEGIN operators
  // BEGIN index operators
  template <class utype> _subtensor<utype> xtensor<utype>::operator()(std::vector<int> indices) const {
#ifdef _CHECK_BOUNDS_
    std::stringstream message;
    uint isize = indices.size();
    if (isize > ndim) {
      message << "Too many indices for tensor with " << ndim << " dimensions (" << isize << " provided).";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    } else {
      for (uint i = 0; i < isize; i++) {
        if ((indices[i] < lindex[i]) || indices[i] > uindex[i]) {
          message << "Index " << indices[i] << " out of bounds for dimension " << i;
          message << " (lindex = " << lindex[i] << ", uindex = " << uindex[i] << ").";
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
      }
    }
#endif
    return _subtensor<utype>(indices, corpus, *this);
  }
#define AST_TEMPLATE(utype) template _subtensor<utype> xtensor<utype>::operator()(std::vector<int>) const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype> xtensor<utype>::operator()(aurostd::xvector<int> indices) const {
#ifdef _CHECK_BOUNDS_
    std::stringstream message;
    if (indices.rows > (int) ndim) {
      message << "Too many indices for tensor with " << ndim << " dimensions (" << indices.rows << " provided).";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    } else {
      for (int i = 0, ixvec = indices.lrows; ixvec <= indices.urows; i++, ixvec++) {
        if (indices[ixvec] < lindex[i] || indices[ixvec] > uindex[i]) {
          message << "Index " << indices[i] << " out of bounds for dimension " << i;
          message << " (lindex = " << lindex[i] << ", uindex = " << uindex[i] << ").";
          throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
      }
    }
#endif
    return _subtensor<utype>(indices, corpus, *this);
  }
#define AST_TEMPLATE(utype) template _subtensor<utype> xtensor<utype>::operator()(aurostd::xvector<int>) const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> _subtensor<utype> xtensor<utype>::operator[](int i) const {
#ifdef _CHECK_BOUNDS_
    if ((i < lindex[0]) || (i > uindex[0])) {
      std::stringstream message;
      message << "Index " << i << " out of bounds for dimension 0";
      message << " (lindex = " << lindex[0] << ", uindex = " << uindex[0] << ").";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
#endif
    return _subtensor<utype>(i, corpus, *this);
  }
#define AST_TEMPLATE(utype) template _subtensor<utype> xtensor<utype>::operator[](int) const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  // END index operators

  //// BEGIN unary operators
  ////// BEGIN with scalars
  template <class utype> xtensor<utype>& xtensor<utype>::operator+=(utype scalar) {
    for (int i = 0; i < nelements; i++) {
      corpus[i] += scalar;
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template xtensor<utype>& xtensor<utype>::operator+=(utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>& xtensor<utype>::operator-=(utype scalar) {
    for (int i = 0; i < nelements; i++) {
      corpus[i] -= scalar;
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template xtensor<utype>& xtensor<utype>::operator-=(utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>& xtensor<utype>::operator*=(utype scalar) {
    for (int i = 0; i < nelements; i++) {
      corpus[i] *= scalar;
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template xtensor<utype>& xtensor<utype>::operator*=(utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>& xtensor<utype>::operator/=(utype scalar) {
    for (int i = 0; i < nelements; i++) {
      corpus[i] /= scalar;
    }
    return *this;
  }
#define AST_TEMPLATE(utype) template xtensor<utype>& xtensor<utype>::operator/=(utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  ////// END with scalars

  ////// BEGIN with xtensors
  template <class utype> xtensor<utype>& xtensor<utype>::operator+=(const xtensor<utype>& tensor) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "operator+= : Dimensions tensor 1:";
      for (uint i = 0; i < ndim; i++) {
        std::cerr << " " << shape[i];
      }
      std::cerr << ", dimensions tensor 2:";
      for (uint i = 0; i < tensor.ndim; i++) {
        std::cerr << " " << tensor.shape[i];
      }
      std::cerr << std::endl;
    }
    if (sameShape(tensor)) {
      for (int i = 0; i < nelements; i++) {
        corpus[i] += tensor.get(i);
      }
      return *this;
        // Exception handling
    } else {
      const std::string message = "Tensors are of different size.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype>& xtensor<utype>::operator+=(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>& xtensor<utype>::operator-=(const xtensor<utype>& tensor) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "operator-= : Dimensions tensor 1:";
      for (uint i = 0; i < ndim; i++) {
        std::cerr << " " << shape[i];
      }
      std::cerr << ", dimensions tensor 2:";
      for (uint i = 0; i < tensor.ndim; i++) {
        std::cerr << " " << tensor.shape[i];
      }
      std::cerr << std::endl;
    }
    if (sameShape(tensor)) {
      for (int i = 0; i < nelements; i++) {
        corpus[i] -= tensor.get(i);
      }
      return *this;
        // Exception handling
    } else {
      const std::string message = "Tensors are of different size.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype>& xtensor<utype>::operator-=(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator+(const xtensor<utype>& tensor) {
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator+(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator-(xtensor<utype> tensor) {
    tensor *= -1;
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator-(const xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  ////// END with xtensors

  ////// BEGIN with _subtensors
  template <class utype> xtensor<utype>& xtensor<utype>::operator+=(const _subtensor<utype>& st) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "operator-= : Dimensions tensor:";
      for (uint i = 0; i < ndim; i++) {
        std::cerr << " " << shape[i];
      }
      std::cerr << ", dimensions subtensor:";
      if (st.indexed_dim == st._tensor.ndim) {
        std::cerr << " 1" << std::endl;
      } else {
        for (uint i = st.indexed_dim; i < st._tensor.ndim; i++) {
          std::cerr << " " << st._tensor.shape[i];
        }
        std::cerr << std::endl;
      }
    }
    if (sameShape(st)) {
      if (st.indexed_dim == st._tensor.ndim) {
        corpus[0] += st.get(st.shift);
      } else {
        const int n = st.getNumElements();
        for (int i = 0, ist = st.shift; i < n; i++, ist++) {
          corpus[i] += st.get(ist);
        }
      }
      return *this;
    } else {
      const std::string message = "Tensor and subtensor are of different size.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype>& xtensor<utype>::operator+=(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype>& xtensor<utype>::operator-=(const _subtensor<utype>& st) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "operator-= : Dimensions tensor:";
      for (uint i = 0; i < ndim; i++) {
        std::cerr << " " << shape[i];
      }
      std::cerr << ", dimensions subtensor:";
      if (st.indexed_dim == st._tensor.ndim) {
        std::cerr << " 1" << std::endl;
      } else {
        for (uint i = st.indexed_dim; i < st._tensor.ndim; i++) {
          std::cerr << " " << st._tensor.shape[i];
        }
        std::cerr << std::endl;
      }
    }
    if (sameShape(st)) {
      if (st.indexed_dim == st._tensor.ndim) {
        corpus[0] -= st.get(st.shift);
      } else {
        const int n = st.getNumElements();
        for (int i = 0, ist = st.shift; i < n; i++, ist++) {
          corpus[i] -= st.get(ist);
        }
      }
      return *this;
    } else {
      const std::string message = "Tensor and subtensor are of different size.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype>& xtensor<utype>::operator-=(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  ////// END with _subtensors
  //// END unary operators

  //// BEGIN binary operators
  ////// BEGIN with scalars
  template <class utype> xtensor<utype> operator*(xtensor<utype> tensor, const utype& scalar) {
    tensor *= scalar;
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator*(xtensor<utype>, const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator*(const utype& scalar, const xtensor<utype>& tensor) {
    return tensor * scalar;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator*(const utype&, const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator/(xtensor<utype> tensor, const utype& scalar) {
    tensor /= scalar;
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator/(xtensor<utype>, const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  ////// END with scalars

  ////// BEGIN with xtensors
  template <class utype> xtensor<utype> operator+(xtensor<utype> tensor1, const xtensor<utype>& tensor2) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "operator+ : Dimensions tensor 1";
      for (uint i = 0; i < tensor1.ndim; i++) {
        std::cerr << " " << tensor1.shape[i];
      }
      std::cerr << ", dimensions tensor 2:";
      for (uint i = 0; i < tensor2.ndim; i++) {
        std::cerr << " " << tensor2.shape[i];
      }
      std::cerr << std::endl;
    }
    if (tensor1.sameShape(tensor2)) {
      tensor1 += tensor2;
      return tensor1;
        // Exception handling
    } else {
      const std::string message = "Tensors are of different size.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator+(xtensor<utype>, const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator-(xtensor<utype> tensor1, const xtensor<utype>& tensor2) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "operator- : Dimensions tensor 1";
      for (uint i = 0; i < tensor1.ndim; i++) {
        std::cerr << " " << tensor1.shape[i];
      }
      std::cerr << ", dimensions tensor 2:";
      for (uint i = 0; i < tensor2.ndim; i++) {
        std::cerr << " " << tensor2.shape[i];
      }
      std::cerr << std::endl;
    }
    if (tensor1.sameShape(tensor2)) {
      tensor1 -= tensor2;
      return tensor1;
        // Exception handling
    } else {
      const std::string message = "Tensors are of different size.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator-(xtensor<utype>, const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  ////// END with xtensors
  ////// BEGIN with _subtensors
  template <class utype> xtensor<utype> operator+(const xtensor<utype>& tensor, const _subtensor<utype>& st) {
    const xtensor<utype> tensor_from_st(st);
    return tensor + tensor_from_st;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator+(const xtensor<utype>&, const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator+(const _subtensor<utype>& st, const xtensor<utype>& tensor) {
    return tensor + st;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator+(const _subtensor<utype>&, const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator+(const _subtensor<utype>& st1, const _subtensor<utype>& st2) {
    const xtensor<utype> tensor1(st1);
    const xtensor<utype> tensor2(st2);
    return tensor1 + tensor2;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator+(const _subtensor<utype>&, const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator-(const xtensor<utype>& tensor, const _subtensor<utype>& st) {
    const xtensor<utype> tensor_from_st(st);
    return tensor - tensor_from_st;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator-(const xtensor<utype>&, const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator-(const _subtensor<utype>& st, const xtensor<utype>& tensor) {
    return -tensor + st;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator-(const _subtensor<utype>&, const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> operator-(const _subtensor<utype>& st1, const _subtensor<utype>& st2) {
    const xtensor<utype> tensor1(st1);
    const xtensor<utype> tensor2(st2);
    return tensor1 - tensor2;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> operator-(const _subtensor<utype>&, const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  ////// END with _subtensors
  //// END binary operators
  // END operators

  // BEGIN member functions
  // set/////////////////////////////////////////////////////////////////////////
  template <class utype> void xtensor<utype>::set(const utype& value) {
    for (int i = 0; i < nelements; i++) {
      corpus[i] = (utype) value;
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::set(const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::set(const utype& value, const int& i) {
    corpus[i] = value;
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::set(const utype&, const int&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::reset() {
    for (int i = 0; i < nelements; i++) {
      corpus[i] = (utype) 0.0;
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::reset();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::clear() {
    reset();
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::clear();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype xtensor<utype>::get(const int& i) const {
    return corpus[i];
  }
#define AST_TEMPLATE(utype) template utype xtensor<utype>::get(const int&) const;
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // END member functions

  // BEGIN non-member functions
  //// BEGIN conversion functions
  ////// From xtensor
  template <class utype> std::vector<utype> xtensor2vector(xtensor<utype>& tensor) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "xtensor2vector - Number of tensor dimensions: " << tensor.ndim << std::endl;
    }
    if (tensor.ndim > 1) {
        // Exception handling
      std::stringstream message;
      message << "Cannot convert xtensor with " << tensor.ndim << " dimensions to vector.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    } else {
      const int rows = tensor.shape[0];
      std::vector<utype> vec(rows);
      for (int i = 0; i < rows; i++) {
        vec[i] = tensor[i + tensor.lindex[0]];
      }
      return vec;
    }
  }
#define AST_TEMPLATE(utype) template std::vector<utype> xtensor2vector(xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> std::vector<utype> xtensor2vector(const _subtensor<utype>& st) {
    xtensor<utype> tensor(st);
    return xtensor2vector(tensor);
  }
#define AST_TEMPLATE(utype) template std::vector<utype> xtensor2vector(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xvector<utype> xtensor2xvector(xtensor<utype>& tensor) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "xtensor2xvector - Number of tensor dimensions: " << tensor.ndim << std::endl;
    }
    if (tensor.ndim > 1) {
        // Exception handling
      std::stringstream message;
      message << "Cannot convert xtensor with " << tensor.ndim << " dimensions to xvector.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    } else {
      const int rows = tensor.shape[0];
      aurostd::xvector<utype> xvec(rows);
      for (int i = 0; i < rows; i++) {
        xvec(i + xvec.lrows) = tensor[i + tensor.lindex[0]];
      }
      return xvec;
    }
  }
#define AST_TEMPLATE(utype) template xvector<utype> xtensor2xvector(xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xvector<utype> xtensor2xvector(const _subtensor<utype>& st) {
    xtensor<utype> tensor(st);
    return xtensor2xvector(tensor);
  }
#define AST_TEMPLATE(utype) template xvector<utype> xtensor2xvector(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xmatrix<utype> xtensor2xmatrix(xtensor<utype>& tensor) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "xtensor2xmatrix - Number of tensor dimensions: " << tensor.ndim << std::endl;
    }
    if (tensor.ndim > 2) {
        // Exception handling
      std::stringstream message;
      message << "Cannot convert xtensor with " << tensor.ndim << " dimensions to xmatrix.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    } else {
      int rows;
      int cols;
      rows = tensor.shape[0];
      if (tensor.ndim > 1) {
        cols = tensor.shape[1];
      } else {
        cols = 1;
      }
      aurostd::xmatrix<utype> mat(rows, cols);
      for (int i = 0; i < rows; i++) {
        if (tensor.ndim > 1) {
          for (int j = 0; j < cols; j++) {
            mat(i + 1, j + 1) = tensor[i + tensor.lindex[0]][j + tensor.lindex[1]];
          }
        } else {
          mat(i, 1) = tensor[i + tensor.lindex[0]];
        }
      }
      return mat;
    }
  }
#define AST_TEMPLATE(utype) template xmatrix<utype> xtensor2xmatrix(xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xmatrix<utype> xtensor2xmatrix(const _subtensor<utype>& st) {
    xtensor<utype> tensor(st);
    return xtensor2xmatrix(tensor);
  }
#define AST_TEMPLATE(utype) template xmatrix<utype> xtensor2xmatrix(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  ////// To xtensor
  template <class utype> xtensor<utype> vector2xtensor(const std::vector<utype>& vec) {
    const std::vector<int> dim(1, vec.size());
    const xtensor<utype> tensor(dim);
    for (size_t i = 0; i < vec.size(); i++) {
      tensor[i + 1] = vec[i];
    }
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> vector2xtensor(const std::vector<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> xvector2xtensor(const aurostd::xvector<utype>& xvec) {
    const std::vector<int> dim(1, xvec.rows);
    const xtensor<utype> tensor(dim);
    for (int i = 0; i < xvec.rows; i++) {
      tensor[i + 1] = xvec(xvec.lrows + i);
    }
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> xvector2xtensor(const aurostd::xvector<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> xmatrix2xtensor(const aurostd::xmatrix<utype>& xmat) {
    std::vector<int> dim(2);
    dim[0] = xmat.rows;
    dim[1] = xmat.cols;
    const xtensor<utype> tensor(dim);
    for (int i = 0; i < xmat.rows; i++) {
      for (int j = 0; j < xmat.cols; j++) {
        tensor[i + 1][j + 1] = xmat(xmat.lrows + i, xmat.lcols + j);
      }
    }
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> xmatrix2xtensor(const aurostd::xmatrix<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  //// END conversion functions

  //// BEGIN return type xtensor
  template <class utype> xtensor<utype> abs(xtensor<utype> tensor) {
    tensor.abs();
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> abs(xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::abs() {
    for (int i = 0; i < nelements; i++) {
      corpus[i] = aurostd::abs(corpus[i]);
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::abs();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> abs(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return abs(tensor_out);
  }
#define AST_TEMPLATE(utype) template xtensor<utype> abs(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> ceil(xtensor<utype> tensor) {
    tensor.ceil();
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> ceil(xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::ceil() {
    for (int i = 0; i < nelements; i++) {
      corpus[i] = std::ceil(corpus[i]);
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::ceil();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> ceil(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return ceil(tensor_out);
  }
#define AST_TEMPLATE(utype) template xtensor<utype> ceil(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> floor(xtensor<utype> tensor) {
    tensor.floor();
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> floor(xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::floor() {
    for (int i = 0; i < nelements; i++) {
      corpus[i] = std::floor(corpus[i]);
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::floor();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> floor(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return floor(tensor_out);
  }
#define AST_TEMPLATE(utype) template xtensor<utype> floor(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> nint(xtensor<utype> tensor) {
    tensor.nint();
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> nint(xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::nint() {
    for (int i = 0; i < nelements; i++) {
      corpus[i] = aurostd::nint(corpus[i]);
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::nint();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> nint(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return nint(tensor_out);
  }
#define AST_TEMPLATE(utype) template xtensor<utype> nint(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> round(xtensor<utype> tensor, const utype& tol) {
    tensor.round(tol);
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> round(xtensor<utype>, const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::round(const utype& tol) {
    for (int i = 0; i < nelements; i++) {
      corpus[i] = aurostd::roundoff(corpus[i], tol);
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::round(const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> round(const _subtensor<utype>& st, const utype& tol) {
    const xtensor<utype> tensor_out(st);
    return round(tensor_out, tol);
  }
#define AST_TEMPLATE(utype) template xtensor<utype> round(const _subtensor<utype>&, const utype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> sign(xtensor<utype> tensor) {
    tensor.sign();
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> sign(xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> void xtensor<utype>::sign() {
    for (int i = 0; i < nelements; i++) {
      corpus[i] = aurostd::sign(corpus[i]);
    }
  }
#define AST_TEMPLATE(utype) template void xtensor<utype>::sign();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xtensor<utype> sign(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return sign(tensor_out);
  }
#define AST_TEMPLATE(utype) template xtensor<utype> sign(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype>
    // Creates an identity tensor of rank n
  xtensor<utype> identity_tensor(const utype& _type, int n) {  // _type is needed to deduce utype
    if (_type) {
      ;
    }  // to suppress compiler warninigs
    xtensor<utype> tensor(n);
    tensor.set((utype) 0.0); // In case the default initialization value is different from 0
    std::vector<int> diag(n);
    for (int i = 1; i <= n; i++) {
      for (int j = 0; j < n; j++) {
        diag[j] = i;
      }
      tensor(diag) = (utype) 1;
    }
    return tensor;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> identity_tensor(const utype&, int);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  //// END return type xtensor

  //// BEGIN return type utype
  template <class utype> utype max(xtensor<utype> tensor) {
    return tensor.max();
  }
#define AST_TEMPLATE(utype) template utype max(xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype xtensor<utype>::max() {
    utype m = corpus[0];
    for (int i = 1; i < nelements; i++) {
      if (corpus[i] > m) {
        m = corpus[i];
      }
    }
    return m;
  }
#define AST_TEMPLATE(utype) template utype xtensor<utype>::max();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype max(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return max(tensor_out);
  }
#define AST_TEMPLATE(utype) template utype max(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype min(xtensor<utype> tensor) {
    return tensor.min();
  }
#define AST_TEMPLATE(utype) template utype min(xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype xtensor<utype>::min() {
    utype m = corpus[0];
    for (int i = 1; i < nelements; i++) {
      if (corpus[i] < m) {
        m = corpus[i];
      }
    }
    return m;
  }
#define AST_TEMPLATE(utype) template utype xtensor<utype>::min();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype min(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return min(tensor_out);
  }
#define AST_TEMPLATE(utype) template utype min(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype sum(xtensor<utype> tensor) {
    return tensor.sum();
  }
#define AST_TEMPLATE(utype) template utype sum(xtensor<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype xtensor<utype>::sum() {
    utype s = 0;
    for (int i = 0; i < nelements; i++) {
      s += corpus[i];
    }
    return s;
  }
#define AST_TEMPLATE(utype) template utype xtensor<utype>::sum();
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype sum(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return sum(tensor_out);
  }
#define AST_TEMPLATE(utype) template utype sum(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype trace(const xtensor<utype>& tensor) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XTENSOR_);
    if (LDEBUG) {
      std::cerr << _XTENSOR_ERR_PREFIX_ << "trace - Tensor shape:";
      for (uint i = 0; i < tensor.ndim; i++) {
        std::cerr << " " << tensor.shape[i];
      }
      std::cerr << std::endl;
    }
    utype trc = 0;
    if (tensor.is_cubic) {
      const uint n = tensor.ndim;
      std::vector<int> diag(n);
      for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++) {
          diag[j] = i + tensor.lindex[j];
        }
        trc += tensor(diag);
      }
      return trc;
    } else {
      const std::string message = "Trace is only defined for cubic tensors.";
      throw xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }
#define AST_TEMPLATE(utype) template utype trace(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype trace(const _subtensor<utype>& st) {
    const xtensor<utype> tensor_out(st);
    return trace(tensor_out);
  }
#define AST_TEMPLATE(utype) template utype trace(const _subtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  //// END return type utype
  // END non-member functions

}  // namespace aurostd
/******************************** END XTENSOR *******************************/

/*********************************** EIJK ***********************************/
// eijk (Ricci Tensor) from old xtensor - used in aflow_xatom.cpp
namespace aurostd {
  int eijk(int i, int j, int k) {
    const int ii = (i - 1) % 3 + 1;
    const int jj = (j - 1) % 3 + 1;
    const int kk = (k - 1) % 3 + 1;
    if ((ii == 1) && (jj == 2) && (kk == 3)) {
      return 1;
    }
    if ((ii == 2) && (jj == 3) && (kk == 1)) {
      return 1;
    }
    if ((ii == 3) && (jj == 1) && (kk == 2)) {
      return 1;
    }
    if ((ii == 1) && (jj == 3) && (kk == 2)) {
      return -1;
    }
    if ((ii == 3) && (jj == 2) && (kk == 1)) {
      return -1;
    }
    if ((ii == 2) && (jj == 1) && (kk == 3)) {
      return -1;
    }
    return 0;
  }

  int eijk(xvector<int> ijk) {
    return eijk(ijk[1], ijk[2], ijk[3]);
  }

  int estarijk(int i, int j, int k) {
    const int ii = (i - 1) % 3 + 1;
    const int jj = (j - 1) % 3 + 1;
    const int kk = (k - 1) % 3 + 1;
    if ((ii == 1) && (jj == 2) && (kk == 3)) {
      return 1;
    }
    if ((ii == 2) && (jj == 3) && (kk == 1)) {
      return 1;
    }
    if ((ii == 3) && (jj == 1) && (kk == 2)) {
      return 1;
    }
    return 0;
  }

  int estarijk(xvector<int> ijk) {
    return estarijk(ijk[1], ijk[2], ijk[3]);
  }
} // namespace aurostd
/*********************************** EIJK ***********************************/

// ----------------------------------------------------------------------------
namespace aurostd {
  /// @brief convert a tensor in Voigt notation to a proper tensor
  ///
  /// @param cij 6x6 tensor in Voigt notation
  ///
  /// @return 3x3x3x3 tensor
  ///
  /// @authors
  /// @mod{SD,20250828,created function}
  template <class utype> xtensor<utype> convert_voigt_to_tensor(const xmatrix<utype>& cij) {
    const xtensor<utype> cijkl(vector<int>{3, 3, 3, 3});
    const std::function<int(int, int)> voigt_imap = [](int i, int j) { return (i == j) ? i : 9 - (i + j); };
    int ii, jj;
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        for (int k = 1; k <= 3; k++) {
          for (int l = 1; l <= 3; l++) {
            ii = voigt_imap(i, j);
            jj = voigt_imap(k, l);
            cijkl[i][j][k][l] = cij[ii][jj];
          }
        }
      }
    }
    return cijkl;
  }
#define AST_TEMPLATE(utype) template xtensor<utype> convert_voigt_to_tensor(const xmatrix<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  /// @brief convert a 3x3x3x3 tensor to a tensor in Voigt notation
  ///
  /// @param 3x3x3x3 tensor
  ///
  /// @return cij 6x6 tensor in Voigt notation
  ///
  /// @authors
  /// @mod{SD,20250828,created function}
  template <class utype> xmatrix<utype> convert_tensor_to_voigt(const xtensor<utype>& cijkl) {
    xmatrix<utype> cij(6, 6);
    vector<std::pair<int, int>> voigt_map = {{1, 1}, {2, 2}, {3, 3}, {2, 3}, {1, 3}, {1, 2}};
    int i, j, k, l;
    for (size_t ii = 1; ii <= 6; ii++) {
      for (size_t jj = 1; jj <= 6; jj++) {
        i = voigt_map[ii - 1].first;
        j = voigt_map[ii - 1].second;
        k = voigt_map[jj - 1].first;
        l = voigt_map[jj - 1].second;
        cij[ii][jj] = cijkl[i][j][k][l];
      }
    }
    return cij;
  }
#define AST_TEMPLATE(utype) template xmatrix<utype> convert_tensor_to_voigt(const xtensor<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
} // namespace aurostd

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************
