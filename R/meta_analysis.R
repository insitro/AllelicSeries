# Purpose: Inverse variance weighted meta-analysis by annotation category.
# Updated: 2024-07-26

#' Inverse Variance Meta-Analysis
#' 
#' @section Notes:
#' This wrapper function is needed because C++ does not allow for default
#' arguments with dynamically allocated values.
#' 
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param beta (snps x 1) vector of effect sizes for 
#'   the coding genetic variants within a gene.
#' @param se (snps x 1) vector of standard errors for the effect sizes.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Defaults to the identity matrix.
#' @param weights Annotation category weights. Defaults to c(1, 2, 3).
#' @return Data.frame with the following columns \itemize{
#' \item anno: gene annotation
#' \item beta_meta: meta-analyzed effect size
#' \item se_meta: standard error of the meta-analyzed effect size
#' }
#' @export
IVWSS <- function(
  anno,
  beta,
  se,
  ld = NULL,
  weights = c(1, 2, 3)
) {
  if (is.null(ld)) {
    n_snp <- length(anno)
    ld <- diag(n_snp)
  }
  out <- IVWCpp(
    anno = anno,
    beta = beta,
    ld = ld,
    se = se,
    weights = weights
  )
  return(out)
}
