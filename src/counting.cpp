// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// For debugging: Rcpp::Rcout << << std::endl; 

// ----------------------------------------------------------------------------

//' Count Variants and Carriers
//'
//' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
//' @param geno (n x snps) genotype matrix.
//' @param min_mac Minimum minor allele count for inclusion. Default: 0.
//' @return Data.frame of allele, variant, and carrier counts.
//' @export
// [[Rcpp::export]]
SEXP Counts(
  arma::colvec anno,
  arma::mat geno,
  const int min_mac = 0
){

  // Apply MAC filter.
  arma::rowvec macs = arma::sum(geno, 0);
  geno = geno.cols(arma::find(macs.t() > min_mac));
  anno = anno.elem(arma::find(macs.t() > min_mac));

  arma::colvec alleles = arma::zeros(3);
  arma::colvec variants = arma::zeros(3);
  arma::colvec carriers = arma::zeros(3);
  
// Loop over annotations.
  for(int j=0; j<3; j++) {
    
    arma::uvec key = arma::find(anno == j);
    arma::mat geno_subset = geno.cols(key);
    variants(j) = geno_subset.n_cols;
    alleles(j) = arma::accu(geno_subset);
    
    // Total number of variants carried by each subject.
    arma::colvec col_sum = arma::sum(geno_subset, 1);
    carriers(j) = arma::accu(col_sum != 0);
    
  };

  return Rcpp::DataFrame::create(
    Rcpp::Named("anno")=arma::linspace(0, 2, 3),
    Rcpp::Named("alleles")=alleles,
    Rcpp::Named("variants")=variants,
    Rcpp::Named("carriers")=carriers
  );
}

