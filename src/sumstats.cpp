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
  const arma::mat x 
){
  return Rcpp::wrap(arma::cor(x));
};


//' Inverse Variance Meta-Analysis
//' 
//' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
//' @param beta (snps x 1) vector of effect sizes for 
//'   the coding genetic variants within a gene.
//' @param ld (snps x snps) matrix of correlations among the genetic variants.
//' @param se (snps x 1) vector of standard errors for the effect sizes.
//' @param weights (3 x 1) vector of annotation category weights.
//' @return Data.frame with the following columns \itemize{
//' \item anno: gene annotation
//' \item beta_meta: meta-analyzed effect size
//' \item se_meta: standard error of the meta-analyzed effect size
//' }
// [[Rcpp::export]]
SEXP IVWCpp(
  const arma::colvec anno,
  const arma::colvec beta,
  const arma::colvec se,
  const arma::mat ld,
  const arma::colvec weights
){

  arma::colvec beta_meta = arma::zeros(3);
  arma::colvec se_meta = arma::zeros(3);
  for(int i=0; i<3; i++){

    // Subset LD matrix.
    arma::uvec key = arma::find(anno == i);
    arma::mat ld_anno = ld.submat(key, key);

    // Calculate variance.
    arma::colvec se_anno = se.elem(key);
    arma::mat d_anno = arma::diagmat(se_anno);
    arma::mat v_anno = d_anno * ld_anno * d_anno;
    arma::mat v_anno_inv = arma::inv_sympd(v_anno);
    double weight_sum = arma::accu(v_anno_inv);

    // Meta-analyze.
    arma::colvec beta_anno = beta.elem(key);
    beta_meta(i) = arma::accu(v_anno_inv * beta_anno) / weight_sum / weights(i);
    se_meta(i) = 1 / sqrt(weight_sum) / weights(i);

  };

  return Rcpp::DataFrame::create(
    Rcpp::Named("anno")=arma::linspace(0, 2, 3),
    Rcpp::Named("beta_meta")=beta_meta,
    Rcpp::Named("se_meta")=se_meta
  );

};


//' Category Correlation
//' 
//' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
//' @param ld (snps x snps) matrix of correlations among the genetic variants.
//' @param maf (snps x 1) vector of minor allele frequencies.
// [[Rcpp::export]]
SEXP CatCor(
  const arma::colvec anno,
  const arma::mat ld,
  const arma::colvec maf
){

  // Initialize to identity matrix.
  arma::mat results = arma::eye(3, 3);
  if (all(maf == 0)) {
    return Rcpp::wrap(results);
  };

  // Only off-diagonals require filling.
  for (int i=0; i<3; i++) {
    for (int j=i+1; j<3; j++) {
      arma::uvec idx1 = arma::find(anno == i);
      arma::uvec idx2 = arma::find(anno == j);

      // Numerator.
      double num = 0;
      for (int m=0; m<idx1.n_elem; m++) {
        for (int n=0; n<idx2.n_elem; n++) {
          double p = maf[idx1[m]];
          double q = maf[idx2[n]];
          num += ld(idx1[m], idx2[n]) * sqrt(4 * p * q * (1 - p) * (1 - q));
        }
      }

      // Denominator term 1.
      double denom1 = 0;
      for (int m=0; m<idx1.n_elem; m++) {
        double p = maf[idx1[m]];
        denom1 += 2 * p * (1 - p);
        for (int n=0; n<idx1.n_elem; n++) {
          if (idx1[n] < idx1[m]) {
            double q = maf[idx1[n]];
            denom1 += 2 * ld(idx1[m], idx1[n]) * sqrt(4 * p * q * (1 - p) * (1 - q));
          }
        }
      }

      // Denominator term 2.
      double denom2 = 0;
      for (int m=0; m<idx2.n_elem; m++) {
        double p = maf[idx2[m]];
        denom2 += 2 * p * (1 - p);
        for (int n=0; n<idx2.n_elem; n++) {
          if (idx2[n] < idx2[m]) {
            double q = maf[idx2[n]];
            denom2 += 2 * ld(idx2[m], idx2[n]) * sqrt(4 * p * q * (1 - p) * (1 - q));
          }
        }
      }

      double score = num / sqrt(denom1 * denom2);
      results(i, j) = results(j, i) = score; 
    }
  }

  return Rcpp::wrap(results);
};


//' Baseline Counts Test from Sumstats
//'
//' @param beta (snps x 1) vector of effect sizes for 
//'   the coding genetic variants within a gene.
//' @param ld (snps x snps) matrix of correlations among the genetic variants.
//' @param se (snps x 1) vector of standard errors for the effect sizes.
//' @return Numeric p-value.
//' @export
// [[Rcpp::export]]
SEXP BaseCountsSS(
  const arma::colvec beta,
  const arma::mat ld,
  const arma::colvec se
){

  // Calculate z-scores.
  const arma::colvec z = beta / se;

  // Inverse variance weights.
  const arma::colvec w = 1 / se;

  // Score statistics.
  const arma::colvec u = w % z;

  // Calculate covariance.
  const arma::mat v = (w * w.t()) % ld;

  // Test statistic.
  arma::mat t = u.t() * arma::solve(v, u, arma::solve_opts::likely_sympd);
  const double tstat = arma::as_scalar(t);

  // P-value.
  int df = u.n_elem;
  Rcpp::Environment base("package:stats");
  Rcpp::Function pval = base["pchisq"]; 
  return pval(Rcpp::_["q"]=tstat, Rcpp::_["df"]=df, Rcpp::_["lower.tail"]=false);
}

