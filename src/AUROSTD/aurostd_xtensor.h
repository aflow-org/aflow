//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************

#ifndef _AUROSTD_XTENSOR_H_
#define _AUROSTD_XTENSOR_H_

namespace aurostd {
  template<class utype> class xtensor; // Forward declaration of xtensor for _subtensor
  template<class utype>
    class _subtensor {
      public:
        _subtensor(const int&, utype*, const xtensor<utype>&);
        _subtensor(const std::vector<int>&, utype*, const xtensor<utype>&);
        _subtensor(const aurostd::xvector<int>&, utype*, const xtensor<utype>&);
        ~_subtensor();

        const xtensor<utype>& _tensor;
        uint indexed_dim;
        int shift;

        _subtensor<utype>& operator[] (const int&);
        _subtensor<utype>& operator() (const std::vector<int>&);
        _subtensor<utype>& operator() (const aurostd::xvector<int>&);
        _subtensor<utype>& operator+= (utype scalar);
        _subtensor<utype>& operator-= (utype scalar);
        _subtensor<utype>& operator*= (utype scalar);
        _subtensor<utype>& operator/= (utype scalar);
        _subtensor<utype>& operator+=(const _subtensor<utype>&);
        _subtensor<utype>& operator-=(const _subtensor<utype>&);
        _subtensor<utype>& operator+=(const xtensor<utype>&);
        _subtensor<utype>& operator-=(const xtensor<utype>&);

        _subtensor<utype>& operator=(const utype&);
        _subtensor<utype>& operator=(const _subtensor<utype>&);
        _subtensor<utype>& operator=(const xtensor<utype>&);
        _subtensor<utype>& operator=(const std::vector<utype>&);
        _subtensor<utype>& operator=(const aurostd::xvector<utype>&);
        _subtensor<utype>& operator=(const aurostd::xmatrix<utype>&);
        operator utype() const;

        int getNumElements() const;
        utype get(const int&) const;
        void set(const utype&);
        void set(const utype&, const int&);

      private:
        bool sameShape(const _subtensor<utype>&);
        bool sameShape(const xtensor<utype>&);
        utype* corpus;
    };
} //namespace aurostd

namespace aurostd {
  template<class utype>
    class xtensor {
      public:
        xtensor();  // Default constructor
        xtensor(const std::vector<int>&);  // Constructor using vectors
        xtensor(const std::vector<int>&, const std::vector<int>&);  // Constructor using vectors
        xtensor(const aurostd::xvector<int>&);  // Constructor using xvectors
        xtensor(const aurostd::xvector<int>&, const aurostd::xvector<int>&);  // Constructor using xvectors
        xtensor(const int&, const int& li = 1);  // Constructor for a tensor with n dimensions of size n
        xtensor(const _subtensor<utype>&); // Constructor from _subtensor
        xtensor(const xtensor<utype>&); // Copy Constructor
        xtensor<utype> operator=(const xtensor<utype>&); // Copy constructor
        xtensor<utype> operator=(const _subtensor<utype>&); // Copy constructor
        ~xtensor();  // Destructor
        void buildTensor(const std::vector<int>&, const std::vector<int>&);

        int* shape; // The shape of the tensor
        bool is_cubic;  // True when all dimensions have the same size
        int nelements;  // Number of elements in the tensor
        uint ndim;  // Number of dimensions of the tensor
        int* uindex;  // The upper indices of each tensor dimension
        int* lindex;  // The lower indices of each tensor dimension
        int* shifts;  // Number of values to skip when incrementing an index in each dimension

        // Indexing operators
        _subtensor<utype> operator()(std::vector<int>) const;
        _subtensor<utype> operator()(aurostd::xvector<int>) const;
        _subtensor<utype> operator[](int) const;

        // Unary operators
        //// with scalars
        xtensor<utype>& operator+=(utype);
        xtensor<utype>& operator-=(utype);
        xtensor<utype>& operator*=(utype);
        xtensor<utype>& operator/=(utype);

        //// with xtensors
        xtensor<utype>& operator+=(const xtensor<utype>&);
        xtensor<utype>& operator-=(const xtensor<utype>&);
        xtensor<utype>& operator+=(const _subtensor<utype>&);
        xtensor<utype>& operator-=(const _subtensor<utype>&);

        // Functions
        void set(const utype&);
        void set(const utype&, const int&);
        void reset();
        void clear();
        bool sameShape(const xtensor<utype>&) const;
        bool sameShape(const _subtensor<utype>&) const;

        utype get(const int&) const;
        void abs();
        void ceil();
        void floor();
        void nint();
        void round(const utype&);
        void sign();
        utype max();
        utype min();
        utype sum();

      private:
        utype *corpus;  // Tensor values
        bool checkInit(const std::vector<int>&, const std::vector<int>&);
        char size;  // Size of the data type used for the tensor
        long int tsize;  // Memory to allocate
    };
} // namespace aurostd

// Unary operators
namespace aurostd {
  template<class utype>
    xtensor<utype> operator+(const xtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator-(xtensor<utype>) __xprototype;
  template<class utype>
    xtensor<utype> operator+(const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator-(const _subtensor<utype>&) __xprototype;
} // namespace aurostd

// Binary operators
namespace aurostd {
  // with scalars
  template<class utype>
    xtensor<utype> operator*(xtensor<utype>, const utype&) __xprototype;
  template<class utype>
    xtensor<utype> operator*(const utype&, const xtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator/(xtensor<utype>, const utype&) __xprototype;
  // with xtensors
  template<class utype>
    xtensor<utype> operator+(xtensor<utype>, const xtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator-(xtensor<utype>, const xtensor<utype>&) __xprototype;
  // with subtensors
  template<class utype>
    xtensor<utype> operator+(const xtensor<utype>&, const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator+(const _subtensor<utype>&, const xtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator+(const _subtensor<utype>&, const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator-(const xtensor<utype>&, const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator-(const _subtensor<utype>&, const xtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> operator-(const _subtensor<utype>&, const _subtensor<utype>&) __xprototype;
} // namespace aurostd

// Functions
namespace aurostd {
  template<class utype>
    std::vector<utype> xtensor2vector(xtensor<utype>&) __xprototype;
  template<class utype>
    std::vector<utype> xtensor2vector(const _subtensor<utype>&) __xprototype;
  template<class utype>
    aurostd::xvector<utype> xtensor2xvector(xtensor<utype>&) __xprototype;
  template<class utype>
    aurostd::xvector<utype> xtensor2xvector(const _subtensor<utype>&) __xprototype;
  template<class utype>
    aurostd::xmatrix<utype> xtensor2xmatrix(xtensor<utype>&) __xprototype;
  template<class utype>
    aurostd::xmatrix<utype> xtensor2xmatrix(const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> vector2xtensor(const std::vector<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> xvector2xtensor(const aurostd::xvector<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> xmatrix2xtensor(const aurostd::xmatrix<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> abs(xtensor<utype>) __xprototype;
  template<class utype>
    xtensor<utype> abs(const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> ceil(xtensor<utype>) __xprototype;
  template<class utype>
    xtensor<utype> ceil(const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> floor(xtensor<utype>) __xprototype;
  template<class utype>
    xtensor<utype> floor(const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> nint(xtensor<utype>) __xprototype;
  template<class utype>
    xtensor<utype> nint(const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> round(xtensor<utype>, const utype& tol=(utype) AUROSTD_ROUNDOFF_TOL) __xprototype;
  template<class utype>
    xtensor<utype> round(const _subtensor<utype>&, const utype& tol=(utype) AUROSTD_ROUNDOFF_TOL) __xprototype;
  template<class utype>
    xtensor<utype> sign(xtensor<utype>) __xprototype;
  template<class utype>
    xtensor<utype> sign(const _subtensor<utype>&) __xprototype;
  template<class utype>
    utype max(xtensor<utype>) __xprototype;
  template<class utype>
    utype max(const _subtensor<utype>&) __xprototype;
  template<class utype>
    utype min(xtensor<utype>) __xprototype;
  template<class utype>
    utype min(const _subtensor<utype>&) __xprototype;
  template<class utype>
    utype sum(xtensor<utype>) __xprototype;
  template<class utype>
    utype sum(const _subtensor<utype>&) __xprototype;
  template<class utype>
    utype trace(const xtensor<utype>&) __xprototype;
  template<class utype>
    utype trace(const _subtensor<utype>&) __xprototype;
  template<class utype>
    xtensor<utype> identity_tensor(const utype&, int) __xprototype;
} // namespace aurostd

// Ricci Tensor functions from old xtensor
namespace aurostd {
  int eijk(int, int, int);
  int eijk(aurostd::xvector<int>);
  int estarijk(int, int, int);
  int estarijk(aurostd::xvector<int>);
}

#endif
//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************
