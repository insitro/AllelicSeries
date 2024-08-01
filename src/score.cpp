// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------

//' Calculate Residual Variance
//'
//' @param y (n x 1) Numeric phenotype vector.
//' @param X (n x q) Numeric covariate matrix.
//'
//' @return Scalar residual variance.
// [[Rcpp::export]]
SEXP ResidVar(
  const arma::colvec y,
  const arma::mat X
){
  const arma::mat YY = y.t() * y;
  const arma::colvec XY = X.t() * y;
  const arma::mat XX = X.t() * X;
  const arma::mat YPY = XY.t() * arma::solve(XX, XY, arma::solve_opts::likely_sympd);
  const double ss = arma::as_scalar(YY - YPY);
  const double out = ss / (X.n_rows - X.n_cols);
  return Rcpp::wrap(out);
}


//' Calculate Score Statistic
//'
//' @param y (n x 1) Numeric phenotype vector.
//' @param G (n x p) Numeric genotype matrix.
//' @param X (n x q) Numeric covariate matrix.
//' @param v Scalar residual variance.
//'
//' @return Scalar score statistic.
// [[Rcpp::export]]
SEXP Score(
  const arma::colvec y,
  const arma::mat G,
  const arma::mat X,
  const double v
){
  
  // Residual.
  const arma::mat XX = X.t() * X;
  const arma::mat XY = X.t() * y;
  const arma::colvec e = y - X * arma::solve(XX, XY, arma::solve_opts::likely_sympd);

  // Score.
  const arma::colvec s = G.t() * e;

  // Calculate G'PG.
  const arma::mat GG = G.t() * G;
  const arma::mat XG = X.t() * G;
  const arma::mat GPG =  GG - XG.t() * arma::solve(XX, XG, arma::solve_opts::likely_sympd);

  // Test statistic.
  const arma::mat stat = s.t() *  arma::solve(GPG, s, arma::solve_opts::likely_sympd);
  const double out = arma::as_scalar(stat / v);

  return Rcpp::wrap(out);
}


