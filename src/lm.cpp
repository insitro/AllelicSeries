// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// ----------------------------------------------------------------------------

//' Ordinary Least Squares
//'
//' Fits the standard OLS model.
//'
//' @param y (n x 1) Numeric vector.
//' @param X (n x p) Numeric matrix.
//'
//' @return List containing the following:
//' \itemize{
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{SE}{Standard errors.}
//' \item{Z}{Z-scores.}
//' }
//' @export
// [[Rcpp::export]]
SEXP OLS(const arma::colvec y, const arma::mat X){
  // Observations
  const int n = y.size();

  // Estimated parameters.
  const int p = X.n_cols;

  // Information.
  const arma::mat A = X.t() * X;

  // Estimate beta.
  const arma::vec b = arma::solve(A, X.t() * y, arma::solve_opts::likely_sympd);

  // Calculate residuals.
  const arma::vec eps = (y - X * b);

  // Scale.
  const double v = arma::as_scalar((eps.t() * eps) / (n - p));

  // Information.
  const arma::mat Ibb = A / v;

  // Standard errors.
  const arma::vec se = arma::sqrt(arma::diagvec(arma::pinv(Ibb)));

  // Z-scores.
  const arma::vec z = b / se;

  return Rcpp::List::create(
    Rcpp::Named("beta") = b,
    Rcpp::Named("v") = v,
    Rcpp::Named("se") = se,
    Rcpp::Named("z") = z,
    Rcpp::Named("resid") = eps
  );
}
