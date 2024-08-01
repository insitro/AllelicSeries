# Purpose: Implement sumstats-based allelic series test.
# Updated: 2024-07-31

# Default weights.
DEFAULT_WEIGHTS <- c(1, 2, 3)

#' Allelic Series Burden Test from Summary Statistics
#' 
#' Allelic series burden test from summary statistics. 
#'
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param beta (snps x 1) vector of effect sizes for the coding genetic variants
#'   within a gene.
#' @param se (snps x 1) vector of standard errors for the effect sizes.
#' @param check Run input checks? Default: TRUE.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Although ideally provided, an identity matrix is assumed if not.
#' @param maf (snps x 1) vector of minor allele frequencies. Although ideally
#'   provided, defaults to the zero vector.
#' @param method Method for aggregating across categories:
#'   {"none", "sum"}. Default: "none".
#' @param weights (3 x 1) vector of annotation category weights.
#' @return Numeric p-value of the allelic series burden test.
#' @examples
#' # Generate data.
#' data <- DGP(n = 1e3)
#' sumstats <- CalcSumstats(data = data)
#' 
#' # Run allelic series burden test from sumstats.
#' results <- ASBTSS(
#'   anno = sumstats$anno,
#'   beta = sumstats$sumstats$beta, 
#'   maf = sumstats$maf,
#'   se = sumstats$sumstats$se,
#'   ld = sumstats$ld
#' )
#' show(results)
#' @export
ASBTSS <- function(
  anno,
  beta,
  se,
  check = TRUE,
  ld = NULL,
  maf = NULL,
  method = "none",
  weights = DEFAULT_WEIGHTS
){
  
  # Check for LD and MAF.
  if (check) {
    CheckInputsSS(
      anno = anno,
      beta = beta,
      se = se,
      maf = maf,
      ld = ld
    )
  }

  n_snps <- length(anno)
  if (is.null(ld)) {ld <- diag(n_snps)}
  if (is.null(maf)) {maf <- rep(0, n_snps)}
  
  # Per-category summary statistics.
  category_sumstats <- IVWSS(
    anno = anno,
    beta = beta,
    se = se,
    ld = ld,
    weights = weights
  )
  
  # Category correlation matrix.
  r3 <- CatCor(anno = anno, ld = ld, maf = maf)
  
  # Run burden test.
  if (method == "none") {
    
    pval <- BaseCountsSS(
      beta = category_sumstats$beta_meta,
      ld = r3,
      se = category_sumstats$se_meta
    )
    
  } else if (method == "sum") {
    
    gene_sumstats <- IVWSS(
      anno = c(0, 0, 0),
      beta = category_sumstats$beta_meta,
      ld = r3,
      se = category_sumstats$se_meta,
      weights = c(1, 1, 1)
    )
    
    pval <- BaseCountsSS(
      beta = gene_sumstats$beta_meta[1],
      se = gene_sumstats$se_meta[1],
      ld = diag(1)
    )
    
  } else {
    
    stop("Select method from among 'none' or 'sum'.")
    
  }
  
  return(pval)
}


# ------------------------------------------------------------------------------

#' Allelic Series SKAT-O from Summary Statistics
#' 
#' Allelic series sequence kernel association test from summary statistics.
#' 
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param beta (snps x 1) vector of effect sizes for the coding genetic variants
#'   within a gene.
#' @param se (snps x 1) vector of standard errors for the effect sizes.
#' @param check Run input checks? Default: TRUE.
#' @param maf (snps x 1) vector of minor allele frequencies. Although ideally
#'   provided, defaults to the zero vector.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Although ideally provided, an identity matrix is assumed if not.
#' @param weights (3 x 1) vector of annotation category weights.
#' @return Numeric p-value of the allelic series SKAT-O test.
#' @examples
#' # Generate data.
#' data <- DGP(n = 1e3)
#' sumstats <- CalcSumstats(data = data)
#' 
#' # Run allelic series SKAT test from sumstats.
#' results <- ASKATSS(
#'   anno = sumstats$anno,
#'   beta = sumstats$sumstats$beta, 
#'   maf = sumstats$maf,
#'   se = sumstats$sumstats$se,
#'   ld = sumstats$ld
#' )
#' show(results)
#' @export
ASKATSS <- function(
    anno, 
    beta, 
    se, 
    check = TRUE,
    ld = NULL,
    maf = NULL,
    weights = DEFAULT_WEIGHTS
){
  
  # Input checks.
  if (check) {
    CheckInputsSS(
      anno = anno,
      beta = beta,
      se = se,
      maf = maf,
      ld = ld
    )
  }
  
  n_snps <- length(anno)
  if (is.null(maf)) {maf <- rep(0, n_snps)}
  
  # Calculate SKAT weights.
  skat_weights <- CalcSkatWeights(
    anno = anno,
    maf = maf,
    weights = weights
  )
  skat_weights <- as.numeric(skat_weights)
  
  # Calculate weighted effect size and standard errors.
  n_snp <- length(beta)
  beta <- beta / skat_weights
  se <- se / skat_weights
  
  # Calculate score test statistic from summary statistics.
  w <- 1 / se
  z <- beta / se
  u <- w * z
  
  # Covariance of Z scores.
  if (!is.null(ld)) {
    cov_z <- (w %*% t(w)) * ld
  } else{
    cov_z <- diag(w^2)
  }
  
  # SKAT-O grid search (provided by SKAT-O paper).
  rhos <- c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 0.999) 
  n_rho <- length(rhos)
  
  # Q statistic: interpolates burden and variance comp statistics.
  qstats <- (1 - rhos) * c(t(u) %*% (u)) + rhos * sum(u)^2
  
  # Loop over values of rho.
  lambdas <- lapply(seq_len(n_rho), function(i) {
    rho <- rhos[i]
    cor_mat <- diag(rep(1 - rho, n_snp)) + matrix(rho, nrow = n_snp, ncol = n_snp)
    cor_mat_sqrt <- chol(cor_mat, pivot = TRUE)
    kernel <- cor_mat_sqrt %*% cov_z %*% t(cor_mat_sqrt)
    out <- as.numeric(GetLambda(kernel))
    return(out)
  })

  # Test parameters.
  opt_params <- SkatOptimalParam(cov_z = cov_z, rhos = rhos)
  results <- PerRhoResults(lambdas = lambdas, qstats = qstats, rhos = rhos)
  pval <- OptimalPval(opt_params = opt_params, results = results)
  return(pval)
}


# ------------------------------------------------------------------------------


#' COding-variant Allelic Series Test from Summary Statistics
#'
#' Main function for performing the allelic series test from summary statistics.
#' Performs both Burden and SKAT type tests, then combines the results to
#' calculate an omnibus p-value. Note that not all tests included in
#' \code{\link{COAST}} are available when working with summary statistics.
#'
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param beta (snps x 1) vector of effect sizes for the coding genetic variants
#'   within a gene.
#' @param se (snps x 1) vector of standard errors for the effect sizes.
#' @param check Run input checks? Default: TRUE.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Although ideally provided, an identity matrix is assumed if not.
#' @param maf (snps x 1) vector of minor allele frequencies. Although ideally
#'   provided, defaults to the zero vector.
#' @param pval_weights (3 x 1) vector of relative weights for combining the
#'   component tests to perform the omnibus test.
#' @param weights (3 x 1) vector of annotation category weights.
#' @return Numeric p-value.
#' @examples 
#' # Generate data.
#' data <- DGP(n = 1e3)
#' sumstats <- CalcSumstats(data = data)
#' 
#' # Run the Coding-variant Allelic Series Test from summary statistics.
#' results <- COASTSS(
#'   anno = sumstats$anno,
#'   beta = sumstats$sumstats$beta, 
#'   maf = sumstats$maf,
#'   se = sumstats$sumstats$se,
#'   ld = sumstats$ld
#' )
#' show(results)
#' @export
COASTSS <- function(
    anno,
    beta, 
    se,
    check = TRUE,
    maf = NULL,
    ld = NULL,
    pval_weights = c(1, 1, 1),
    weights = DEFAULT_WEIGHTS
) {
  
  # Input checks.
  if (check) {
    CheckInputsSS(
      anno = anno,
      beta = beta,
      se = se,
      maf = maf,
      ld = ld
    )
  }
  
  # Baseline model p-value.
  p_base <- ASBTSS(
    anno = anno,
    beta = beta,
    se = se,
    check = FALSE,
    ld = ld,
    maf = maf,
    method = "none",
    weights = c(1, 1, 1)
  )
    
  # Sum model p-value.
  p_burden <- ASBTSS(
    anno = anno,
    beta = beta,
    se = se,
    check = FALSE,
    ld = ld,
    maf = maf,
    method = "sum",
    weights = weights
  )
  
  # SKAT model p-value.
  p_skat <- ASKATSS(
    anno = anno,
    beta = beta,
    se = se,
    check = FALSE,
    ld = ld,
    maf = maf,
    weights = weights
  )

  # Omnibus p-value.
  pvals <- c(p_base, p_burden, p_skat)
  p_omni <- RNOmni::OmniP(pvals, pval_weights)
  
  # Output.
  df_pvals <- data.frame(
    test = c("baseline", "sum_count", "allelic_skat", "omni"),
    type = c("burden", "burden", "skat", "omni"),
    pval = c(pvals, p_omni)
  )
  
  # Output.
  out <- methods::new(
    Class = "COAST",
    Counts = NULL,
    Pvals = df_pvals
  )
  return(out)
}
