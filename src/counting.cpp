// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------

//' Count Variants and Carriers
//'
//' @param anno (snps x 1) annotation vector with integer values in 1 through
//'   the number of annotation categories L.
//' @param geno (n x snps) genotype matrix.
//' @param n_anno Number of annotation categories L.
//' @param min_mac Minimum minor allele count for inclusion. Default: 0.
//' @return Data.frame of allele, variant, and carrier counts.
//' @export
// [[Rcpp::export]]
SEXP Counts(
  arma::colvec anno,
  arma::mat geno,
  const int n_anno,
  const int min_mac = 0
){

  // Apply MAC filter.
  arma::rowvec macs = arma::sum(geno, 0);
  geno = geno.cols(arma::find(macs.t() > min_mac));
  anno = anno.elem(arma::find(macs.t() > min_mac));

  arma::colvec alleles = arma::zeros(n_anno);
  arma::colvec variants = arma::zeros(n_anno);
  arma::colvec carriers = arma::zeros(n_anno);
  
  // Loop over annotations.
  for (int j=0; j<n_anno; j++) {
    
    arma::uvec key = arma::find(anno == (j + 1));
    arma::mat geno_subset = geno.cols(key);
    variants(j) = geno_subset.n_cols;
    alleles(j) = arma::accu(geno_subset);
    
    // Total number of variants carried by each subject.
    arma::colvec col_sum = arma::sum(geno_subset, 1);
    carriers(j) = arma::accu(col_sum != 0);
    
  };

  return Rcpp::DataFrame::create(
    Rcpp::Named("anno")=arma::linspace(1, n_anno, n_anno),
    Rcpp::Named("alleles")=alleles,
    Rcpp::Named("variants")=variants,
    Rcpp::Named("carriers")=carriers
  );
}

