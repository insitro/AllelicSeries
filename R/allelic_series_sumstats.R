# Purpose: Implement sumstats-based allelic series test.
# Updated: 2024-11-13

#' Allelic Series Burden Test from Summary Statistics
#' 
#' Allelic series burden test from summary statistics. 
#' 
#' @section Notes:
#' * The allelic series burden does not require the minor allele frequencies.
#'
#' @param anno (snps x 1) annotation vector with integer values in 1 through
#'   the number of annotation categories L.
#' @param beta (snps x 1) vector of effect sizes for the coding genetic variants
#'   within a gene.
#' @param se (snps x 1) vector of standard errors for the effect sizes.
#' @param check Run input checks? Default: TRUE.
#' @param eps Epsilon added to the diagonal of the LD matrix if not positive
#'   definite. Note, smaller values increase the chances of a false positive. 
#' @param lambda Optional genomic inflation factor. Defaults to 1, which
#'   results in no rescaling.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Although ideally provided, an identity matrix is assumed if not.
#' @param method Method for aggregating across categories:
#'   ("none", "sum"). Default: "none".
#' @param weights (L x 1) vector of annotation category weights. Note that the
#'   number of annotation categories L is inferred from the length of `weights`.
#' @return Numeric p-value of the allelic series burden test.
#' @examples
#' # Generate data.
#' data <- DGP(n = 1e3)
#' sumstats <- CalcSumstats(data = data)
#' 
#' # Run allelic series burden test from sumstats.
#' results <- ASBTSS(
#'   anno = sumstats$sumstats$anno,
#'   beta = sumstats$sumstats$beta, 
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
  eps = 1,
  lambda = 1,
  ld = NULL,
  method = "none",
  weights = c(1, 2, 3)
){
  
  # Ensure annotations are one-indexed.
  anno <- RelevelAnno(anno)
  
  # Check for LD and MAF.
  if (check) {
    is_pd <- CheckInputsSS(
      anno = anno,
      beta = beta,
      se = se,
      lambda = lambda,
      ld = ld,
      weights = weights,
      is_skat = FALSE
    )
  } else if (!is.null(ld)) {
    is_pd <- isPD(ld)
  } else{
    is_pd <- TRUE
  }
  
  n_snps <- length(anno)
  n_anno <- length(weights)
  if (is.null(ld)) {ld <- diag(n_snps)}
  if (!is_pd) {ld <- ld + eps * diag(n_snps)}
  
  # Run burden test.
  if (method == "none") {
    
    pval <- BaselineSS(
      anno = anno,
      beta = beta,
      ld = ld,
      se = se,
      n_anno = n_anno
    )
    
  } else if (method == "sum") {
    
    pval <- SumCountSS(
      anno = anno,
      beta = beta,
      ld = ld,
      se = se,
      weights = weights
    )
    
  } else {
    
    stop("Select method from among 'none' or 'sum'.")
    
  }
  
  # Apply genomic control.
  lambda <- max(1, lambda)
  pval <- GenomicControl(lambda = lambda, pval = pval)
  return(pval)
}


# ------------------------------------------------------------------------------

#' Allelic Series SKAT-O from Summary Statistics
#' 
#' Allelic series sequence kernel association test from summary statistics.
#' 
#' @section Notes: 
#' * The SKAT test requires per-variant minor allele frequencies (MAFs) for 
#'   the purpose of up-weighting rarer variants. If unknown, `maf` can be 
#'   safely omitted. The only consequence is that the SKAT weights will no 
#'   longer be inversely proportional to the genotypic variance.
#' 
#' @param anno (snps x 1) annotation vector with integer values in 1 through
#'   the number of annotation categories L.
#' @param beta (snps x 1) vector of effect sizes for the coding genetic variants
#'   within a gene.
#' @param se (snps x 1) vector of standard errors for the effect sizes.
#' @param check Run input checks? Default: TRUE.
#' @param eps Epsilon added to the diagonal of the LD matrix if not positive
#'   definite. Note, smaller values increase the chances of a false positive. 
#' @param lambda Optional genomic inflation factor. Defaults to 1, which
#'   results in no rescaling.
#' @param maf (snps x 1) vector of minor allele frequencies. Although ideally
#'   provided, defaults to the zero vector.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Although ideally provided, an identity matrix is assumed if not.
#' @param weights (L x 1) vector of annotation category weights. Note that the
#'   number of annotation categories L is inferred from the length of `weights`.
#' @return Numeric p-value of the allelic series SKAT-O test.
#' @examples
#' # Generate data.
#' data <- DGP(n = 1e3)
#' sumstats <- CalcSumstats(data = data)
#' 
#' # Run allelic series SKAT test from sumstats.
#' # Note: the SKAT test requires MAF.
#' results <- ASKATSS(
#'   anno = sumstats$sumstats$anno,
#'   beta = sumstats$sumstats$beta, 
#'   maf = sumstats$sumstats$maf,
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
    eps = 1,
    lambda = 1,
    ld = NULL,
    maf = NULL,
    weights = c(1, 2, 3)
){
  
  # Ensure annotations are one-indexed.
  anno <- RelevelAnno(anno)
  
  # Input checks.
  if (check) {
    sink <- CheckInputsSS(
      anno = anno,
      beta = beta,
      se = se,
      lambda = lambda,
      ld = ld,
      weights = weights,
      is_skat = TRUE,
      maf = maf
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
  n_snps <- length(beta)
  beta <- beta / skat_weights
  se <- se / skat_weights
  
  # Calculate score test statistic from summary statistics.
  w <- 1 / se
  z <- beta / se
  u <- w * z
  
  # Covariance of Z scores.
  if (!is.null(ld)) {
    cov_z <- (w %*% t(w)) * ld
    if (!isPD(cov_z)) {cov_z <- cov_z + eps * diag(n_snps)}
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
    cor_mat <- diag(rep(1 - rho, n_snps)) + 
      matrix(rho, nrow = n_snps, ncol = n_snps)
    cor_mat_sqrt <- chol(cor_mat, pivot = TRUE)
    kernel <- cor_mat_sqrt %*% cov_z %*% t(cor_mat_sqrt)
    out <- as.numeric(GetLambda(kernel))
    return(out)
  })

  # Test parameters.
  opt_params <- SkatOptimalParam(cov_z = cov_z, rhos = rhos)
  results <- PerRhoResults(lambdas = lambdas, qstats = qstats, rhos = rhos)
  pval <- OptimalPval(opt_params = opt_params, results = results)
  
  # Apply genomic control.
  lambda <- max(1, lambda)
  pval <- GenomicControl(lambda = lambda, pval = pval)
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
#' @param anno (snps x 1) annotation vector with integer values in 1 through
#'   the number of annotation categories L.
#' @param beta (snps x 1) vector of effect sizes for the coding genetic variants
#'   within a gene.
#' @param se (snps x 1) vector of standard errors for the effect sizes.
#' @param check Run input checks? Default: TRUE.
#' @param eps Epsilon added to the diagonal of the LD matrix if not positive
#'   definite. Note, epsilon should increase as the sample size decreases. 
#' @param lambda Optional (3 x 1) vector of inflation factors, one for each
#'   component test. Defaults to a 1s vector, which results in no rescaling.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Although ideally provided, an identity matrix is assumed if not.
#' @param maf (snps x 1) vector of minor allele frequencies. Although ideally
#'   provided, defaults to the zero vector.
#' @param pval_weights (3 x 1) vector of relative weights for combining the
#'   component tests to perform the omnibus test. The default of
#'   c(0.25, 0.25, 0.50) gives the SKAT test equal weight to the two burden tests. 
#' @param weights (L x 1) vector of annotation category weights. Note that the
#'   number of annotation categories L is inferred from the length of `weights`.
#' @return Numeric p-value.
#' @examples 
#' # Generate data.
#' data <- DGP(n = 1e3)
#' sumstats <- CalcSumstats(data = data)
#' 
#' # Run the Coding-variant Allelic Series Test from summary statistics.
#' results <- COASTSS(
#'   anno = sumstats$sumstats$anno,
#'   beta = sumstats$sumstats$beta, 
#'   se = sumstats$sumstats$se,
#'   ld = sumstats$ld,
#'   maf = sumstats$sumstats$maf,
#' )
#' show(results)
#' @export
COASTSS <- function(
    anno,
    beta, 
    se,
    check = TRUE,
    eps = 1,
    lambda = c(1, 1, 1),
    ld = NULL,
    maf = NULL,
    pval_weights = c(0.25, 0.25, 0.50),
    weights = c(1, 2, 3)
) {
  
  # Ensure annotations are one-indexed.
  anno <- RelevelAnno(anno)
  
  # Input checks.
  if (check) {
    sink <- CheckInputsSS(
      anno = anno,
      beta = beta,
      se = se,
      lambda = lambda,
      ld = ld,
      weights = weights,
      maf = maf,
      is_skat = TRUE
    )
  }
  
  # Baseline model p-value.
  n_anno <- length(weights)
  unif_weights <- rep(1, n_anno)
  p_base <- ASBTSS(
    anno = anno,
    beta = beta,
    se = se,
    check = FALSE,
    eps = eps,
    ld = ld,
    method = "none",
    weights = unif_weights
  )
    
  # Sum model p-value.
  p_burden <- ASBTSS(
    anno = anno,
    beta = beta,
    se = se,
    check = FALSE,
    eps = eps,
    ld = ld,
    method = "sum",
    weights = weights
  )
  
  # SKAT model p-value.
  p_skat <- tryCatch({
    ASKATSS(
      anno = anno,
      beta = beta,
      se = se,
      check = FALSE,
      eps = eps,
      ld = ld,
      maf = maf,
      weights = weights
    )
  },
    error = function(cond) {return(NA)}
  )

  # Genomic control.
  lambda <- pmax(lambda, 1)
  p_base <- GenomicControl(lambda[1], p_base)
  p_burden <- GenomicControl(lambda[2], p_burden)
  p_skat <- GenomicControl(lambda[3], p_skat)
  
  # Omnibus p-value.
  pvals <- c(p_base, p_burden, p_skat)
  key <- sapply(pvals, function(x) {!is.na(x)})
  p_omni <- RNOmni::OmniP(pvals[key], pval_weights[key])
  
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
