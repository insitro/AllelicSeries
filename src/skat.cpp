// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------

//' Calculate SKAT Weights
//' @param anno (snps x 1) annotation vector with integer values in 1 through
//'   the number of annotation categories L.
//' @param maf (snps x 1) vector of minor allele frequencies.
//' @param weights (L x 1) vector of annotation category weights. Note that the
//'   number of annotation categories L is inferred from the length of `weights`.
//' @return (snps x 1) vector of weights for SKAT test.
//' @noRd 
// [[Rcpp::export]]
SEXP CalcSkatWeights(
  const arma::vec& anno, 
  const arma::vec& maf, 
  const arma::vec& weights
){
  const int n = anno.n_elem;
  arma::colvec skat_weights = arma::zeros(n);
  
  const int n_anno = weights.n_elem;
  for (int j=0; j<n_anno; j++) {
    arma::uvec ids = arma::find(anno == (j + 1));
    skat_weights.elem(ids).fill(weights[j]);
  }
  
  arma::vec v = maf % (1 - maf);
  v.elem(arma::find(v == 0)).fill(1);
  skat_weights = arma::sqrt(skat_weights / v);
  
  return Rcpp::wrap(skat_weights);
}


// ----------------------------------------------------------------------------

//' Get SKAT Eigenvalues
//' @param kernel Symmetric matrix.
//' @return Vector of eigenvalues.
//' @noRd
// [[Rcpp::export]]
SEXP GetLambda(const arma::mat& kernel) {
  // Calculate eigenvalues.
  arma::colvec eigenvalues;
  arma::eig_sym(eigenvalues, kernel);
  
  // Filter eigenvalues (following SKAT paper).
  arma::uvec idx1 = arma::find(eigenvalues >= 0); 
  double mean_threshold = arma::mean(eigenvalues(idx1)) / 100000; 
  arma::uvec idx2 = arma::find(eigenvalues > mean_threshold); 
  arma::colvec out = eigenvalues.elem(idx2);
  return Rcpp::wrap(out);
}


// ----------------------------------------------------------------------------

//' SKAT Optimal Parameters
//' @param cov_z Covariance of the Z-scores
//' @param rhos Vector of rho values to evaluate.
//' @return List of optimal parameter values.
//' @noRd 
// [[Rcpp::export]]
SEXP SkatOptimalParam(const arma::mat& cov_z, const arma::vec& rhos){
  int p_m = cov_z.n_cols;
  int r_n = rhos.n_elem;

  // Compute z_mean_sq
  double z_mean_sq = arma::accu(cov_z) / (p_m * p_m);

  // Compute cof1
  arma::vec cof1 = arma::sum(cov_z, 0).t() / (p_m * z_mean_sq);

  // Compute W3.2.t
  arma::mat w3_2_t = cov_z - cof1 * cof1.t() * z_mean_sq;

  // Compute eigenvalues from W3.2.t
  arma::vec eigenvalues;
  arma::eig_sym(eigenvalues, w3_2_t);

  // Compute W3.3.item (VarRemain)
  double var_remain = arma::accu(cof1 * cof1.t() * z_mean_sq % w3_2_t) * 4;

  // Compute mu_q, var_q, ker_q, and df
  double mu_q = arma::accu(eigenvalues);
  double var_q = arma::accu(arma::square(eigenvalues)) * 2 + var_remain;
  double ker_q = (12 * arma::accu(arma::pow(eigenvalues, 4))) / std::pow(arma::accu(arma::square(eigenvalues)), 2);
  double df = 12 / ker_q;

  // Compute tau for each r.corr in r_all.
  arma::colvec tau = arma::zeros(r_n);
  for (int i = 0; i < r_n; i++) {
      double r_corr = rhos(i);
      double term1 = p_m * p_m * r_corr + arma::accu(arma::square(cof1)) * (1 - r_corr);
      tau(i) = term1 * z_mean_sq;
  }

  // Return the result as List
  return Rcpp::List::create(
    Rcpp::Named("mu_q") = mu_q,
    Rcpp::Named("var_q") = var_q,
    Rcpp::Named("ker_q") = ker_q,
    Rcpp::Named("lambda") = eigenvalues,
    Rcpp::Named("var_remain") = var_remain,
    Rcpp::Named("df") = df,
    Rcpp::Named("tau") = tau
  );
};
