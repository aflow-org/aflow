// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XFIT_CPP_
#define _AUROSTD_XFIT_CPP_

#ifndef _AUROSTD_XFIT_H_
#include "aurostd_xfit.h"
#endif

#define DEBUG_XFIT false // toggles debug output for xfit

//********************************************************************************
// Functions to work with polynomials
//********************************************************************************
namespace aurostd {
  /// @brief evaluates the value of the polynomial with coefficients p at the value x
  ///
  /// @param x value at which to evaluates the polynomial
  /// @param p coefficients of the polynomial
  ///
  /// @return polynomial evaluated at x
  ///
  /// @authors
  /// @mod{AS,20200825,created function}
  /// @mod{SD,20220718,created documentation}
  ///
  /// @note The polynomial is represented in the following form:
  /// p(x) = p0+x*(p1+x*(p2+x*(p3+...x*(p(n-1)+x*(pn)))))
  template<class utype> utype evalPolynomial(const utype x, const xvector<utype>& p) {
    utype res = p[p.urows];
    for (int i=p.urows-1; i>=p.lrows; i--) res = res*x + p[i];
    return res;
  }

  /// @brief evaluates the value of the polynomial with coefficients p at the values x
  ///
  /// @param x values at which to evaluates the polynomial
  /// @param p coefficients of the polynomial
  ///
  /// @return polynomial evaluated at x
  ///
  /// @authors
  /// @mod{SD,20220422,created function}
  template<class utype> xvector<utype> evalPolynomial_xv(const xvector<utype>& x, const xvector<utype>& p) {
    xvector<utype> res(x.rows);
    for (int i = res.lrows; i <= res.urows; i++) {
      res(i) = evalPolynomial(x(i), p);
    }
    return res;
  }

  /// @brief evaluates the value of the polynomial with coefficients p at the values x
  ///
  /// @param x values at which to evaluate the polynomial
  /// @param p coefficients of the polynomial
  ///
  /// @return polynomial evaluated at x
  ///
  /// @authors
  /// @mod{SD,20220505,created function}
  template<class utype> xmatrix<utype> evalPolynomial_xm(const xmatrix<utype>& x, const xvector<utype>& p) {
    xmatrix<utype> res(x.rows, x.cols);
    for (int i = res.lrows; i <= res.urows; i++) {
      for (int j = res.lcols; j <= res.ucols; j++) {
        res(i, j) = evalPolynomial(x(i, j), p);
      }
    }
    return res;
  }

  /// @brief evaluates the value and the derivatives of the polynomial with coefficients p at
  /// the value x and outputs the derivatives in the dp array
  ///
  /// @param x values at which to evaluate the polynomial
  /// @param p coefficients of the polynomial
  /// @param dp coefficients of the derivative of the polynomial multiplied by the monomial base
  ///
  /// @authors
  /// @mod{AS,20200825,created function}
  /// @mod{SD,20220718,created documentation}
  ///
  /// @note The highest derivative is determined by the size of the dp array and the result
  /// is stored in ascending order starting from the zeroth derivative (function itself).
  template<class utype> void evalPolynomialDeriv(const utype x, const xvector<utype>& p, xvector<utype>& dp) {
    for (int i=dp.lrows; i<=dp.urows; i++) dp[i] = 0.0;
    dp[dp.lrows] = p[p.urows];
    for (int i=p.urows-1; i>=p.lrows; i--){
      for (int j=dp.urows; j>=dp.lrows+1; j--) dp[j] = dp[j]*x + dp[j-1];
      dp[dp.lrows] = dp[dp.lrows]*x + p[i];
    }

    int factor = 1;
    for (int i=dp.lrows+1; i<dp.urows; i++){
      factor *= i;
      dp[i+1] *= factor;
    }
  }

  /// @brief evaluates the value and the derivatives of the polynomial with coefficients p at
  /// the value x and returns the result as an array (dp)
  ///
  /// @param x values at which to evaluate the polynomial
  /// @param p coefficients of the polynomial
  /// @param n order of the derivative
  ///
  /// @return Coefficients of the nth derivative of the polynomial multiplied by the monomial base.
  ///
  /// @authors
  /// @mod{AS,20200825,created function}
  /// @mod{SD,20220718,created documentation}
  ///
  /// @note The result is stored in ascending order starting from the zero's derivative
  /// (function itself).
  template<class utype> xvector<utype> evalPolynomialDeriv(const utype x, const xvector<utype>& p, const uint n) {
    xvector<utype> dp(n+1);
    evalPolynomialDeriv(x, p, dp);
    return dp;
  }

  /// @brief returns the coefficients of the nth derivative of the polynomial
  ///
  /// @param x values at which to evaluate the polynomial
  /// @param p coefficients of the polynomial
  /// @param n order of the derivative
  ///
  /// @return coefficients of the nth derivative of the polynomial.
  ///
  /// @authors
  /// @mod{SD,SD20220422,created function}
  template<class utype> xvector<utype> evalPolynomialCoeff(const xvector<utype>& p, const uint n) {
    if ((int)n >= p.urows) {return xvector<utype>(1);}
    xvector<utype> dp(p.urows - n);
    for (int i = p.urows; i > (int)n; i--) {
      dp(i - n) = p(i) * aurostd::factorial<utype>(i - 1) / aurostd::factorial<utype>(i - n - 1);
    }
    return dp;
  }

  /// @brief Constructs the Vandermonde matrix with columns up to a given number of columns n,
  /// which corresponds to the polynomial of degree n-1.
  ///
  /// @param x Values at which to evaluate the polynomial.
  /// @param n Order of the polynomial.
  ///
  /// @return Vandermonde matrix.
  ///
  /// @authors
  /// @mod{AS,20200825,created function}
  /// @mod{SD,20220718,created documentation}
  ///
  /// Vandermonde matrix:
  ///  1 x1^1 x1^2 .. x1^(n-1)
  ///  1 x2^1 x2^2 .. x2^(n-1)
  ///  ...
  ///  1 xm^1 xm^2 .. xm^(n-1)
  /// where m is the size of x array
  template<class utype> xmatrix<utype> Vandermonde_matrix(const xvector<utype>& x, const int n) {
    xmatrix<utype> VM(x.rows, n);

    for (int i=x.lrows; i<=x.urows; i++){
      VM[i][1] = 1.0;
      for (int j=2; j<=n; j++) VM[i][j] = std::pow(x[i], j-1);
    }

    return VM;
  }

  /// @brief calculates the extremum of the polynomial bounded by the region [xmin,xmax] by searching
  /// for the value x where the derivative of polynomial equals to zero using the bisection method
  ///
  /// @param p coefficients of the polynomial
  /// @param xmin minimum value of range to search for the extremum
  /// @param xmax maximum value of range to search for the extremum
  /// @param tol tolerence when to terminate the search
  ///
  /// @return value at which the polynomial has an extremum
  ///
  /// @authors
  /// @mod{AS,20200825,created function}
  /// @mod{SD,20220427,modified function to throw error rather than return zero when convergence cannot be reached}
  ///
  /// @note It is assumed that only one extremum exists between xmin and xmax.
  /// The specific use of this function is to search for the equilibrium volume of a given
  /// polynomial E(V) dependency.
  ///
  /// @note The main idea: if a continuous function f(x) has the root f(x0)=0 in the interval [a,b], 
  /// then its sign in the interval [a,x0) is opposite to its sign in the interval (x0,b].
  /// The root searching algorithm is iterative: one bisects the interval [a,b] into two 
  /// subintervals and for the next iteration picks a subinterval where f(x) has 
  /// opposite sign on its ends.
  /// This process continues until the length of the interval is less than some negligibly
  /// small number and one picks the middle of the interval as x0.
  /// In this situation, f(x0)~0 and this condition is used as a criterion to stop the
  /// iteration procedure.
  template<class utype> utype polynomialFindExtremum(const xvector<utype>& p, const utype xmin, const utype xmax, const utype tol) {
    bool LDEBUG = (FALSE || DEBUG_XFIT || XHOST.DEBUG);
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " begin" << std::endl;

    utype left_end = xmin; utype right_end = xmax;
    utype middle = 0.5*(left_end + right_end);

    xvector<utype> dp(2); // index 1 - value, index 2 - first derivative
    evalPolynomialDeriv(left_end, p, dp);
    utype f_at_left_end  = dp[2];
    evalPolynomialDeriv(right_end, p, dp);
    utype f_at_right_end = dp[2];
    evalPolynomialDeriv(middle, p, dp);
    utype f_at_middle    = dp[2];

    // no root within a given interval
    if (sign(f_at_left_end) == sign(f_at_right_end)) {
      string message = "No root within the given interval";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    if (LDEBUG){
      cerr << __AFLOW_FUNC__ << " left_end= "  << left_end  << "middle= ";
      cerr << middle  << "right_end= "  << right_end  << std::endl;
      cerr << __AFLOW_FUNC__ << " f_left= " << f_at_left_end;
      cerr << "f_middle= " << f_at_middle << "f_right= ";
      cerr << f_at_right_end << std::endl;
    }

    // Iterate until the convergence criterion is reached:
    // f(middle) is sufficiently close to zero.
    // Meanwhile, do a sanity check that the function has opposite signs at the interval ends.
    while ((sign(f_at_left_end) != sign(f_at_right_end)) && (aurostd::abs(f_at_middle)>tol)){
      if (sign(f_at_left_end) == sign(f_at_middle)){
        std::swap(left_end, middle);
        std::swap(f_at_left_end, f_at_middle);
      }

      if (sign(f_at_right_end) == sign(f_at_middle)){
        std::swap(right_end, middle);
        std::swap(f_at_right_end, f_at_middle);
      }

      middle = 0.5*(left_end + right_end);
      evalPolynomialDeriv(middle, p, dp);
      f_at_middle = dp[2];

      if (LDEBUG){
        cerr << __AFLOW_FUNC__ << " left_end= "  << left_end  << "middle= ";
        cerr << middle  << "right_end= "  << right_end  << std::endl;
        cerr << __AFLOW_FUNC__ << " f_left= " << f_at_left_end;
        cerr << "f_middle= " << f_at_middle << "f_right= ";
        cerr << f_at_right_end << std::endl;
      }
    }

    // utype-check that the convergence criterion was reached
    if (aurostd::abs(f_at_middle) > tol) {
      string message = "Bisection method did not converge";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    if (LDEBUG) cerr << __AFLOW_FUNC__ << " end" << std::endl;
    return middle;
  }

  /// @brief calculates the coefficients for a polynomial p(x) of degree n that fits the y data with weights w
  ///
  /// @param x input value of the data points
  /// @param y output value of the data points
  /// @param n order of the polynomial
  /// @param w weight of the data points
  /// @param scale_input scale the input data for the fit
  ///
  /// @return coefficients of the polynomial
  ///
  /// @authors
  /// @mod{SD,20220422,created function}
  ///
  /// @note The weight w_i determines how much one wants to weight the point y_i when performing the fit
  template<class utype> xvector<utype> polynomialCurveFit(const xvector<utype>& _x, const xvector<utype>& _y, const int n, const xvector<utype>& _w, const bool scale_input) {
    bool LDEBUG = (FALSE || DEBUG_XFIT || XHOST.DEBUG);
    xvector<utype> x = scale_input ? (_x - aurostd::mean(_x)) / aurostd::stddev(_x) : _x;
    utype wtot = 0.0;
    for (int i = _w.lrows; i <= _w.urows; i++) {
      if (_w(i) < (utype)0.0) {
        stringstream message;
        message << "Negative weight i=" << i;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      else if (std::isnan(_w(i))) {
        stringstream message;
        message << "NaN weight i=" << i;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      else if (aurostd::abs(_w(i)) == (utype)INFINITY) {
        stringstream message;
        message << "Inf weight i=" << i;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      wtot += _w(i);
    }
    if (aurostd::isequal(wtot, (utype)0.0)) {
      string message = "All the weights cannot be zero";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    xmatrix<utype> VM = Vandermonde_matrix(x, n + 1);
    if (LDEBUG) {cerr << "VM old=" << VM << endl;}
    xvector<utype> w = _w / wtot; // normalize weights
    for (int i = VM.lrows; i <= VM.urows; i++) {
      VM.setmat(w(i) * VM.getxmat(i, i, 1, VM.cols), i, 1);
    }
    if (LDEBUG) {cerr << "VM new=" << VM << endl;}
    xvector<utype> y = aurostd::elementwise_product(w, _y);
    return aurostd::LinearLeastSquares(VM, y);
  }

  /// @brief calculates the (transpose) n-by-n companion matrix of a polynomial, where n is the polynomial degree
  /// transposed companion matrix
  ///
  /// @param p coefficients of the polynomial
  ///
  /// @return companion matrix of the polynomial
  ///
  /// @authors
  /// @mod{SD,20220318,created function}
  ///
  ///   0   1   0  ..  0
  ///   0   0   1  ..  0
  ///  ...
  ///   0   0   0  ..  1
  ///  -p0 -p1 -p2 .. -p(n-1)
  template<class utype> xmatrix<utype> companion_matrix(const xvector<utype>& p) {
    if (p(p.urows) == 0.0) {
      string message = "Leading polynomial coefficient is zero";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    xmatrix<utype> CM(p.urows - 1, p.urows - 1, p.lrows, p.lrows);
    for (int i = p.lrows; i < p.urows - 1; i++) {CM(i, i + 1) = 1.0;}
    for (int i = p.lrows; i <= p.urows - 1; i++) {CM(p.urows - 1, i) = -p(i) / p(p.urows);}
    return CM;
  }

  /// @brief calculates the roots of a polynomial as the eigenvalues of the companion matrix
  ///
  /// @param p coefficients of the polynomial
  /// @param rr real part of the polynomial roots
  /// @param ri imaginary part of the polynomial roots
  ///
  /// @authors
  /// @mod{SD,20220318,created function}
  ///
  /// @note The real and imaginary parts of the roots are returned in rr and ri, respectively.
  ///
  /// @see
  /// @doi{10.1090/S0025-5718-1995-1262279-2}
  template<class utype> void polynomialFindRoots(const xvector<utype>& p, xvector<utype>& rr, xvector<utype>& ri) {
    for (int i = p.lrows; i <= p.urows; i++) {
      if (std::isnan(p(i))) {
        stringstream message;
        message << "NaN in polynomial coefficient i=" << i;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
      else if (aurostd::abs(p(i)) == (utype)INFINITY) {
        stringstream message;
        message << "Inf in polynomial coefficient i=" << i;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
      }
    }
    aurostd::eigen(companion_matrix(p), rr, ri);
  }
}

//********************************************************************************

//********************************************************************************
// Root-finding algorithms
//********************************************************************************
namespace aurostd {
  /// @brief finds the zeros of a univariate function by Brent's method
  ///
  /// @param a lower (or upper) bound of the search space
  /// @param b upper (or lower) bound of the search space
  /// @param f univeriate function, stored as a lambda expression
  /// @param zero value of the argument at which the function is zero
  /// @param niter number of iterations to perform
  /// @param tol tolerance when to terminate the search
  ///
  /// @return boolean whether the zero was found
  ///
  /// @authors
  /// @mod{SD,20220517,created function}
  ///
  /// @see
  /// @xlink{https://people.math.sc.edu/Burkardt/cpp_src/brent/brent.html}
  bool findZeroBrent(const double a, const double b, const std::function<double(double)>& f, double& zero, const uint niter, const double _tol) {
    double c, d, e, fa, fb, fc, m, p, q, r, s, sa, tol;
    if (aurostd::isequal(aurostd::sign(f(a)), aurostd::sign(f(b)))) {
      stringstream message;
      message << "Function evaluation at the endpoints have the same sign, F(a)=" << f(a) << " F(b)=" << f(b);
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    sa = a;
    zero = b;
    fa = f(sa);
    fb = f(zero);
    c = sa;
    fc = fa;
    e = zero - sa;
    d = e;
    uint iter = 0;
    while (iter < niter) {
      if (std::abs(fc) < std::abs(fb)) {
        sa = zero;
        zero = c;
        c = sa;
        fa = fb;
        fb = fc;
        fc = fa;
      }
      tol = 2.0 * _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ * std::abs(zero) + _tol;
      m = 0.5 * (c - zero);
      if (std::abs(m) <= tol || aurostd::isequal(fb, 0.0)) {break;}
      if (std::abs(e) < tol || std::abs(fa) <= std::abs(fb)) {
        e = m;
        d = e;
      }
      else {
        s = fb / fa;
        if (aurostd::isequal(sa, c)) {
          p = 2.0 * m * s;
          q = 1.0 - s;
        }
        else {
          q = fa / fc;
          r = fb / fc;
          p = s * (2.0 * m * q * (q - r) - (zero - sa) * (r - 1.0));
          q = (q - 1.0) * (r - 1.0) * (s - 1.0);
        }
        if (0.0 < p) {
          q = - q;
        }
        else {
          p = - p;
        }
        s = e;
        e = d;
        if (2.0 * p < 3.0 * m * q - std::abs(tol * q) && p < std::abs(0.5 * s * q)) {
          d = p / q;
        }
        else {
          e = m;
          d = e;
        }
      }
      sa = zero;
      fa = fb;
      if (tol < std::abs(d)) {
        zero = zero + d;
      }
      else if (0.0 < m) {
        zero = zero + tol;
      }
      else {
        zero = zero - tol;
      }
      fb = f(zero);
      if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)) {
        c = sa;
        fc = fa;
        e = zero - sa;
        d = e;
      }
      iter++;
    }
    return (iter == niter) ? false : true;
  }

  /// @brief finds a set of zeros for a square system (N variables, N equations) using the Newton-Raphson method
  ///
  /// @param x0 initial guess for the zeros
  /// @param vf vector of functions, stored as lambda expressions
  /// @param jac vector of vector of the first derivatives of the functions (Jacobian), stored as lambda expressions
  /// @param x set of zeros of the system
  ///
  /// @return boolean whether a set zeros was found
  ///
  /// @authors
  /// @mod{SD,20220616,created function}
  ///
  /// @see
  /// @xlink{Appendix A in Power System Operations,https://link.springer.com/book/10.1007/978-3-319-69407-8}
  bool findZeroNewtonRaphson(const xvector<double>& _x0, const vector<std::function<double(xvector<double>)>>& vf, const vector<vector<std::function<double(xvector<double>)>>>& jac, xvector<double>& x, const uint niter, const double tol) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    uint n = (uint)_x0.rows;
    if (n != vf.size() || n != jac[0].size() || n != jac.size()) {
      string message;
      message = "System of equations is not square";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    xvector<double> x0 = _x0, f(n);
    xmatrix<double> J(n, n), iJ, I = aurostd::identity(J);
    bool converged = false;
    uint iter = 0;
    while (iter < niter) {
      for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++) {
          J(i + 1, j + 1) = jac[i][j](x0);
        }
      }
      try { // check for singular matrix
        iJ = aurostd::inverse(J);
      }
      catch (aurostd::xerror& err) {
        J = J + tol * I; // Tikhanov regularization
        iJ = aurostd::inverse(J);
      }
      for (uint i = 0; i < n; i++) {f(i + 1) = vf[i](x0);}
      x = x0 - iJ * f;
      converged = true;
      for (int i = x.lrows; i <= x.urows && converged; i++) {
        if (std::abs(x(i) - x0(i)) > tol) {converged = false;}
      }
      if (converged) {break;}
      x0 = x;
      iter++;
    }
    if (LDEBUG) {cerr << __AFLOW_FUNC__ << " iter=" << iter << " | dx=" << aurostd::modulus(x - x0) << endl;} 
    return (iter == niter) ? false : true;
  }

  /// @brief finds the set of zeros for a square system (N variables, N equations) using gradient deflation
  ///
  /// @param x0 initial guess for the zeros
  /// @param vf vector of functions, stored as lambda expressions
  /// @param jac vector of vector of the first derivatives of the functions (Jacobian), stored as lambda expressions
  /// @param mx matrix of zeros for the system, each column is a different set of zeros
  ///
  /// @return boolean whether a set of zeros was found
  ///
  /// @authors
  /// @mod{SD,20220619,created function}
  ///
  /// @see
  /// @doi{10.1007/bf02165004}
  bool findZeroDeflation(const xvector<double>& x0, const vector<std::function<double(xvector<double>)>>& _vf, const vector<vector<std::function<double(xvector<double>)>>>& _jac, xmatrix<double>& mx, const uint niter, const double tol) {
    mx.clear();
    vector<vector<double>> vx;
    vector<std::function<double(xvector<double>)>> vf = _vf, vdf;
    vector<vector<std::function<double(xvector<double>)>>> jac = _jac;
    vector<std::function<double(xvector<double>)>> gfac;
    vector<vector<std::function<double(xvector<double>)>>> vgfac;
    vector<vector<std::function<double(xvector<double>)>>> dgfac;
    vector<vector<vector<std::function<double(xvector<double>)>>>> vdgfac;
    std::function<double(uint i, xvector<double>)> coeff;
    std::function<double(uint i, uint j, xvector<double>)> dcoeff;
    xvector<double> r;
    while (aurostd::findZeroNewtonRaphson(x0, vf, jac, r, niter, tol)) {
      vx.push_back(aurostd::xvector2vector<double>(r));
      vf.clear();
      vdf.clear();
      jac.clear();
      gfac.clear();
      for (uint i = 0; i < _vf.size(); i++) {
        gfac.push_back(
                        [i, r, _jac](xvector<double> x) {
                          double gfac = 0.0; 
                          for (uint j = 0; j < _jac.size(); j++) {
                            gfac += _jac[i][j](r) * (x(j + 1) - r(j + 1));
                          }
                          return std::pow(gfac, -1.0);
                        }
                      );
      }
      vgfac.push_back(gfac);
      coeff = [vgfac](uint i, xvector<double> x) {
                double coeff = 1.0;
                for (uint ix = 0; ix < vgfac.size(); ix++) {coeff *= vgfac[ix][i](x);} 
                return coeff;
              };
      dgfac.clear();
      for (uint i = 0; i < _vf.size(); i++) {
        for (uint j = 0; j < _vf.size(); j++) {
          vdf.push_back(
                         [i, j, r, _jac, gfac](xvector<double> x) {
                           return -1.0 * _jac[i][j](r) * std::pow(gfac[i](x), 2.0);
                         }
                       );
        }
        dgfac.push_back(vdf);
        vdf.clear();
      }
      vdgfac.push_back(dgfac);
      dcoeff = [coeff, vgfac, vdgfac](uint i, uint j, xvector<double> x) {
                 double dcoeff = 0.0;
                 for (uint ix = 0; ix < vdgfac.size(); ix++) {dcoeff += vdgfac[ix][i][j](x) / vgfac[ix][i](x);} 
                 return coeff(i, x) * dcoeff;
               };
      for (uint i = 0; i < _vf.size(); i++) {
        vf.push_back([i, _vf, coeff](xvector<double> x) {return coeff(i, x) * _vf[i](x);});
        for (uint j = 0; j < _vf.size(); j++) {
          vdf.push_back(
                         [i, j, _vf, _jac, coeff, dcoeff](xvector<double> x) {
                           return coeff(i, x) * _jac[i][j](x) + dcoeff(i, j, x)  * _vf[i](x);
                         }
                       );
        }
        jac.push_back(vdf);
        vdf.clear();
      }
    }
    if (vx.size() == 0) {return false;}
    mx = aurostd::trasp(aurostd::vectorvector2xmatrix<double>(vx));
    return true;
  }
}

//********************************************************************************
// Auxiliary functions
//********************************************************************************
namespace aurostd {
  /// @brief checks whether the numerical and analytical derivatives match
  ///
  /// @param x test value for the arguments at which to valuate the function and its derivatives
  /// @param vf function of interest, stored as a lambda expression
  /// @param df vector of the first derivative of the function, stored as lambda expressions
  /// @param tol tolerance at which the numerical and analytical derivatives are equal
  ///
  /// @return boolean whether the numerical and analytical derivatives are equal
  ///
  /// @authors
  /// @mod{SD,20220619,created function}
  bool checkDerivatives(const xvector<double>& _x, const std::function<double(xvector<double>)>& f, const vector<std::function<double(xvector<double>)>>& df, const double tol) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if ((uint)_x.rows != df.size()) {
      string message = "Number of variables and derivatives must be equal";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    double dx = 10.0 * tol, df_approx, df_exact;
    xvector<double> x1 = _x, x2 = _x;
    bool match = true;
    for (int i = _x.lrows; i <= _x.urows; i++) {
      x1(i) -= dx;
      x2(i) += dx;
      df_approx = 0.5 * (f(x2) - f(x1)) / dx;
      df_exact = df[(uint)i - (uint)_x.lrows](_x);
      if (LDEBUG) {cerr << __AFLOW_FUNC__ << " i=" << i << " | df_approx=" << df_approx << " | df_exact=" << df_exact << endl;} 
      if (!aurostd::isequal(df_approx, df_exact, tol)) {match = false;}
      x1 = _x;
      x2 = _x;
    }
    return match;
  }

  /// @brief numerically calculates the Jacobian of a system of functions using finite difference
  ///
  /// @param vf vector of functions, stored as lambda expressions
  /// @param dx step sizes to perform the central difference
  ///
  /// @return vector of vector of the first derivatives of the functions (Jacobian), stored as lambda expressions
  ///
  /// @authors
  /// @mod{SD,20220624,created function}
  vector<vector<std::function<double(xvector<double>)>>> calcNumericalJacobian(const vector<std::function<double(xvector<double>)>>& vf, const xvector<double>& _dx) {
    vector<std::function<double(xvector<double>)>> df;
    vector<vector<std::function<double(xvector<double>)>>> jac;
    xvector<double> dx = aurostd::abs(_dx);
    for (uint i = 0; i < vf.size(); i++) {
      for (int j = dx.lrows; j <= dx.urows; j++) {
        df.push_back(
                      [i, j, dx, vf](xvector<double> x) {
                        xvector<double> x1 = x, x2 = x; x1(j) -= dx(j); x2(j) += dx(j); 
                        return 0.5 * (vf[i](x2) - vf[i](x1)) / dx(j);
                      }
                    );
      }
      jac.push_back(df);
      df.clear();
    }
    return jac;
  }
}

//********************************************************************************
//              Definitions of the NonlinearFit class members
namespace aurostd{
  NonlinearFit::NonlinearFit()
  {
    free();
  }

  NonlinearFit::NonlinearFit(const NonlinearFit &nlf)
  {
    if (this==&nlf) return;
    free(); copy(nlf);
  }

  NonlinearFit::NonlinearFit(xvector<double>& x_in, xvector<double>& y_in,
      xvector<double>& guess_in,
      double foo(const double x, const xvector<double>& p, xvector<double>& dydp),
      double tol_in, double tau_in, int max_iter_in)
  {
    free();
    Npoints = x_in.rows;
    Nparams = guess_in.rows;
    tol = tol_in;
    tau = tau_in;
    max_iter = max_iter_in;
    x = x_in;
    y = y_in;
    residuals = xvector<double> (Npoints);
    guess = guess_in;
    p = xvector<double> (Nparams);
    dydp = xvector<double> (Nparams);
    A = xmatrix<double> (Nparams, Nparams);
    J = xmatrix<double> (Npoints, Nparams);
    f = foo;
  }

  NonlinearFit::~NonlinearFit() { free(); }

  const NonlinearFit& NonlinearFit::operator=(const NonlinearFit &nlf)
  {
    copy(nlf);
    return *this;
  }

  void NonlinearFit::clear() { free(); }
  void NonlinearFit::free()
  {
    Npoints = 0;
    Nparams = 0;
    tol = 0.0;
    tau = 0.0;
    max_iter = 0;
    x.clear();
    y.clear();
    residuals.clear();
    guess.clear();
    p.clear();
    dydp.clear();
    A.clear();
    J.clear();
    f = NULL;
  }

  void NonlinearFit::copy(const NonlinearFit &nlf)
  {
    if (this==&nlf) return;

    Npoints = nlf.Npoints;
    Nparams = nlf.Nparams;
    tol = nlf.tol;
    tau = nlf.tau;
    max_iter = nlf.max_iter;
    x = nlf.x;
    y = nlf.y;
    residuals = nlf.residuals;
    guess = nlf.guess;
    p = nlf.p;
    dydp = nlf.dydp;
    A = nlf.A;
    J = nlf.J;
    f = nlf.f;
  }

  /// Calculates the residual sum of squares of a model function w.r.t the given data
  ///
  double NonlinearFit::calculateResidualSquareSum(const xvector<double> &params)
  {
    double chi_sqr = 0.0;

    calculateResiduals(params); // calculate residuals and save them to residuals
    for (int i=1; i<= Npoints; i++) chi_sqr += std::pow(residuals[i], 2);

    return chi_sqr;
  }

  /// Calculates residuals of a given model function and stores the results
  ///
  void NonlinearFit::calculateResiduals(const xvector<double> &params)
  {
    for (int i=1; i<=Npoints; i++) residuals[i] = y[i] - f(x[i], params, dydp);
  }

  /// Calculates the Jacobian of a given model function for the given "guess" parameters
  ///
  void NonlinearFit::Jacobian(const xvector<double> &guess)
  {
    for (int i=1; i<=Npoints; i++){
      f(x[i], guess, dydp); 
      for (int j=1; j<=Nparams; j++){
        J[i][j] = -dydp[j]; // derivative of a given function w.r.t fit parameters
      }
    }
  }

  /// Nonlinear fit using the Levenberg-Marquardt algorithm.
  /// The implementation here is based on the ideas from Numerical Recipes and
  /// K.Madsen et al. Methods For Non-linear Least Squares Problems
  /// http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
  ///
  bool NonlinearFit::fitLevenbergMarquardt()
  {
    bool LDEBUG = (FALSE || DEBUG_XFIT || XHOST.DEBUG);
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " begin"  << std::endl;

    static const double step_scaling_factor = 10.0;

    xvector<double> pnew(Nparams);
    xmatrix<double> G(Nparams, 1), M(Nparams, Nparams);

    int iter = 0;
    double chi_sqr = 0.0, new_chi_sqr =0.0; // old and new residual square sums

    p = guess; Jacobian(p); calculateResiduals(p); // p are the fitted parameters

    xmatrix<double> A = trasp(J)*J;
    xvector<double> g = -trasp(J)*residuals;

    // determine initial step
    double lambda = tau*maxDiagonalElement(A);
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " lambda = " << lambda << std::endl;

    while (iter<max_iter){
      iter++;
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " iteration: " << iter << std::endl;

      // M = A + lambda*diag(A)
      M = A; for (int i=1; i<=Nparams; i++) M[i][i] += lambda*A[i][i];

      if (LDEBUG) cerr << __AFLOW_FUNC__ << " M:" << std::endl << M << std::endl;

      // transformation: xvector(N) => xmatrix(N,1) to use GaussJordan function
      for (int i=1; i<=Nparams; i++) G[i][1] = g[i];
      aurostd::GaussJordan(M,G);

      pnew = p + M*g;

      chi_sqr = calculateResidualSquareSum(p);
      new_chi_sqr = calculateResidualSquareSum(pnew);

      if (std::abs(new_chi_sqr - chi_sqr) < tol){
        p = pnew;
        break;
      }

      // update only if step leads to a smaller residual sum of squares
      if (new_chi_sqr < chi_sqr){
        p = pnew; Jacobian(p); calculateResiduals(p);

        A = trasp(J)*J;
        g = -trasp(J)*residuals;

        lambda /= step_scaling_factor;
      }
      else{
        lambda *= step_scaling_factor;
      }
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " pnew = " << p << std::endl;
    }

    if (LDEBUG) cerr << __AFLOW_FUNC__ << " end"  << std::endl;
    return iter < max_iter;
  }
}

#endif
// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
