// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XFIT_H_
#define _AUROSTD_XFIT_H_

//********************************************************************************
// Functions to work with polynomials
//********************************************************************************
namespace aurostd {
  template<class utype> utype evalPolynomial(const utype x, const xvector<utype>& p) __xprototype;
  template<class utype> xvector<utype> evalPolynomial_xv(const xvector<utype>& x, const xvector<utype>& p) __xprototype;
  template<class utype> xmatrix<utype> evalPolynomial_xm(const xmatrix<utype>& x, const xvector<utype>& p) __xprototype;
  template<class utype> void evalPolynomialDeriv(const utype x, const xvector<utype>& p, xvector<utype>& dp) __xprototype;
  template<class utype> xvector<utype> evalPolynomialDeriv(const utype x, const xvector<utype>& p, const uint n) __xprototype;
  template<class utype> xvector<utype> evalPolynomialCoeff(const xvector<utype>& p, const uint n) __xprototype; //SD20220425
  template<class utype> xmatrix<utype> Vandermonde_matrix(const xvector<utype>& x, const int n) __xprototype;
  template<class utype> utype polynomialFindExtremum(const xvector<utype>& p, const utype xmin, const utype xmax,
      const utype tol=(utype)AUROSTD_IDENTITY_TOL) __xprototype;  //CO20230313 - utype clang warnings
  template<class utype> xvector<utype> polynomialCurveFit(const xvector<utype>& x, const xvector<utype>& y, const int n, const xvector<utype>& w, const bool scale_input=false) __xprototype; //SD20220422
  template<class utype> xmatrix<utype> companion_matrix(const xvector<utype>& p) __xprototype; //SD20220318
  template<class utype> void polynomialFindRoots(const xvector<utype>& p, xvector<utype>& rr, xvector<utype>& ri) __xprototype; //SD20220318
}

//********************************************************************************
// Root-finding algorithms
//********************************************************************************
namespace aurostd {
  bool findZeroBrent(const double a, const double b, const std::function<double(double)>& f, double& zero, const uint niter=100, const double tol=(double)AUROSTD_IDENTITY_TOL); //SD20220517 //CO20230313 - utype clang warnings
  bool findZeroNewtonRaphson(const xvector<double>& x0, const vector<std::function<double(xvector<double>)>>& vf, const vector<vector<std::function<double(xvector<double>)>>>& jac, xvector<double>& x, const uint niter=100, const double tol=(double)AUROSTD_IDENTITY_TOL); //SD20220616 //CO20230313 - utype clang warnings
  bool findZeroDeflation(const xvector<double>& x0, const vector<std::function<double(xvector<double>)>>& vf, const vector<vector<std::function<double(xvector<double>)>>>& jac, xmatrix<double>& mx, const uint niter=100, const double tol=(double)AUROSTD_IDENTITY_TOL); //SD20220619  //CO20230313 - utype clang warnings
}

//********************************************************************************
// Auxiliary functions
//********************************************************************************
namespace aurostd {
  bool checkDerivatives(const xvector<double>& x, const std::function<double(xvector<double>)>& f, const vector<std::function<double(xvector<double>)>>& df, const double tol=(double)AUROSTD_IDENTITY_TOL); //SD20220622 //CO20230313 - utype clang warnings
  vector<vector<std::function<double(xvector<double>)>>> calcNumericalJacobian(const vector<std::function<double(xvector<double>)>>& vf, const xvector<double>& _dx); //SD20220624
}

//********************************************************************************
//              Definitions of the NonlinearFit class members
namespace aurostd{
  /// Fit to a nonlinear model using the Levenberg-Marquardt algorithm.
  ///
  /// The implementation here is based on the ideas from Numerical Recipes and
  /// K. Madsen et al. Methods For Non-linear Least Squares Problems
  /// http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
  ///
  /// Caution: the default value for the parameter tau (a scaling factor for the initial step size)
  /// was picked to yield a correct fit to the Murnaghan equation of state. 
  /// If you observe that this is not a good choice for your function,
  /// try 1e-6 if the initial guess is believed to be a good approximation to the true
  /// parameters. Otherwise 1e-3 or even 1 might be a better choice.
  class NonlinearFit{
    public:
      NonlinearFit();
      NonlinearFit(const NonlinearFit &nlf);
      NonlinearFit(xvector<double>& x, xvector<double>& y, xvector<double>& guess,
          double foo(const double x, const xvector<double>& p, xvector<double>& dydp),
          double tol=1e-6, double tau=1e-12, int max_iter=1000);
      ~NonlinearFit();
      const NonlinearFit& operator=(const NonlinearFit &qha);
      int Npoints, Nparams;
      double tol; /// convergence tolerance criterion
      double tau; /// scaling parameter for initial step size
      int max_iter; /// maximum number of allowed iterations
      xvector<double> x,y;   // data points
      xvector<double> residuals; // residuals of a given model function
      xvector<double> guess; // initial guess for fit parameters
      xvector<double> p;     // parameters obtained by fit
      xvector<double> dydp;  // derivative of a given function w.r.t parameters
      xmatrix<double> A;     // J^T.J matrix
      xmatrix<double> J;     // Jacobian of a model function w.r.t parameters
      double (*f)(const double x, const xvector<double> &p, xvector<double> &dydp);
      bool fitLevenbergMarquardt();
      void Jacobian(const xvector<double> &guess);
      void calculateResiduals(const xvector<double> &params);
      double calculateResidualSquareSum(const xvector<double> &params);
      void clear();

    private:
      void free();
      void copy(const NonlinearFit &nlf);
  };
}

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

