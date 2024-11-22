// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath> 
#include <RcppArmadillo.h>
#include <Rcpp.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------

//' Correlation C++
//' 
//' @section Notes:
//' Verified this function is faster that R's built-in correlation function
//' for large genotype matrices.
//'
//' @param x Numeric matrix.
//' @return Numeric matrix of correlation among the columns.
// [[Rcpp::export]]
SEXP CorCpp(
  const arma::mat &x 
){
  return Rcpp::wrap(arma::cor(x));
};


//' Check if Positive Definite
//'
//' @param x Numeric matrix.
//' @param force_symmetry Force the matrix to be symmetric?
//' @param tau Threshold the minimum eigenvalue must exceed for the matrix
//'   to be considered positive definite. 
//' @return Logical indicating whether the matrix is PD.
//' @export 
// [[Rcpp::export]]
SEXP isPD(
  const arma::mat &x, 
  bool force_symmetry=false,
  double tau=1e-8
) {
  
  // Symmetrize.
  arma::mat y = x;
  if (force_symmetry) {
    y = 0.5 * (x + x.t());
  } 

  // Check eigenvalues.
  bool out = false;
  if (y.is_symmetric()) {
      arma::colvec lambda;
      arma::eig_sym(lambda, y);
      out = arma::all(lambda > tau);
  }

  return Rcpp::wrap(out);
};


// Annotation matrix
// 
// @param anno (snps x 1) annotation vector.
// @param n_anno Number of annotation categories L.
// [[Rcpp::export]]
arma::mat AnnoMat(
  const arma::colvec &anno,
  const int n_anno = 3
) {

  arma::mat out = arma::zeros(anno.n_elem, n_anno);
  for (int j=0; j<anno.n_elem; j++) {
    out(j, anno(j) - 1) = 1.0;
  };

  return out;
};


//' Baseline Counts Test from Sumstats
//'
//' @param anno (snps x 1) annotation vector with integer values in 1 through
//'   the number of annotation categories L.
//' @param beta (snps x 1) vector of effect sizes for 
//'   the coding genetic variants within a gene.
//' @param ld (snps x snps) matrix of correlations among the genetic variants.
//' @param se (snps x 1) vector of standard errors for the effect sizes.
//' @param n_anno Number of annotation categories L.
//' @param return_beta Return estimated effect sizes and standard errors?
//'   Default: FALSE.
//' @return If `return_beta`, a list containing the category effect sizes,
//'   standard errors, and the p-value. Otherwise, the numeric p-value only.
//' @export
// [[Rcpp::export]]
SEXP BaselineSS(
  const arma::colvec &anno,
  const arma::colvec &beta,
  const arma::mat &ld,
  const arma::colvec &se,
  const int n_anno = 3,
  const bool return_beta = false
){

  // Filter to require non-zero standard error.
  arma::uvec key = arma::find(se > 0);
  if (key.empty()) {return Rcpp::wrap(0.0);}
  const arma::colvec beta_nz = beta.elem(key);
  const arma::mat ld_nz = ld.submat(key, key);
  const arma::colvec se_nz = se.elem(key);

  // Calculate score statistics.
  const arma::colvec v = se_nz % se_nz;
  const arma::colvec eta = beta_nz / v;
  
  // Category score statistics.
  const arma::mat d = AnnoMat(anno, n_anno);
  const arma::mat u = d.t() * eta;
  
  // Calculate covariance matrix.
  const arma::mat inv_se = arma::diagmat(1 / se);
  const arma::mat cov = d.t() * inv_se * ld * inv_se * d;

  // Test statistic.
  const arma::mat inv_cov = arma::pinv(cov);
  arma::mat t = u.t() * inv_cov * u;
  const double tstat = arma::as_scalar(t);

  // P-value.
  int df = u.n_elem;
  Rcpp::Environment base("package:stats");
  Rcpp::Function pchisq = base["pchisq"]; 
  SEXP pval = pchisq(Rcpp::_["q"]=tstat, Rcpp::_["df"]=df, Rcpp::_["lower.tail"]=false);

  if (!return_beta) {return pval;};

  // Calculate betas and standard errors.
  const arma::colvec beta_hat = inv_cov * u;
  const arma::colvec beta_hat_se = arma::sqrt(arma::diagvec(inv_cov));
  
  // Output.
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_hat,
    Rcpp::Named("se") = beta_hat_se,
    Rcpp::Named("pval") = pval
  );
};


//' Allelic Sum Test from Sumstats
//'
//' @param anno (snps x 1) annotation vector with integer values in 1 through
//'   the number of annotation categories L.
//' @param beta (snps x 1) vector of effect sizes for 
//'   the coding genetic variants within a gene.
//' @param ld (snps x snps) matrix of correlations among the genetic variants.
//' @param se (snps x 1) vector of standard errors for the effect sizes.
//' @param weights (L x 1) vector of annotation category weights. Note that the
//'   number of annotation categories L is inferred from the length of `weights`.
//' @param return_beta Return estimated effect sizes and standard errors?
//'   Default: FALSE.
//' @return If `return_beta`, a list containing the category effect sizes,
//'   standard errors, and the p-value. Otherwise, the numeric p-value only.
//' @export
// [[Rcpp::export]]
SEXP SumCountSS(
  const arma::colvec &anno,
  const arma::colvec &beta,
  const arma::mat &ld,
  const arma::colvec &se,
  const arma::colvec &weights,
  const bool return_beta = false
){

  // Alias.
  const arma::colvec w = weights;
  const int n_anno = weights.n_elem;
  
  // Filter to require non-zero standard error.
  arma::uvec key = arma::find(se > 0);
  if (key.empty()) {return Rcpp::wrap(0.0);}
  const arma::colvec beta_nz = beta.elem(key);
  const arma::mat ld_nz = ld.submat(key, key);
  const arma::colvec se_nz = se.elem(key);

  // Calculate score statistics.
  const arma::colvec v = se_nz % se_nz;
  const arma::colvec eta = beta_nz / v;
  
  // Category score statistics.
  const arma::mat d = AnnoMat(anno, n_anno);
  const arma::mat u = w.t() * d.t() * eta;
  
  // Calculate covariance matrix.
  const arma::mat inv_se = arma::diagmat(1 / se);
  const arma::mat cov = w.t() * d.t() * inv_se * ld * inv_se * d * w;

  // Test statistic.
  const arma::mat inv_cov = arma::pinv(cov);
  arma::mat t = u.t() * inv_cov * u;
  const double tstat = arma::as_scalar(t);

  // P-value.
  int df = u.n_elem;
  Rcpp::Environment base("package:stats");
  Rcpp::Function pchisq = base["pchisq"]; 
  SEXP pval = pchisq(Rcpp::_["q"]=tstat, Rcpp::_["df"]=df, Rcpp::_["lower.tail"]=false);

  if (!return_beta) {return pval;};

  // Calculate betas and standard errors.
  const arma::colvec beta_hat = inv_cov * u;
  const arma::colvec beta_hat_se = arma::sqrt(arma::diagvec(inv_cov));
  
  // Output.
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_hat,
    Rcpp::Named("se") = beta_hat_se,
    Rcpp::Named("pval") = pval
  );

};
