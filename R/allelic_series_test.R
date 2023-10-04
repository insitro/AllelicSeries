# Purpose: Allelic series test.
# Updated: 2023-10-04

# Default weights.
DEFAULT_WEIGHTS <- c(1, 2, 3)

#' Aggregator
#'
#' Aggregates genotypes within annotation categories.
#'
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param geno (n x snps) genotype matrix.
#' @param drop_empty Drop empty columns? Default: TRUE.
#' @param indicator Convert raw counts to indicators? Default: FALSE.
#' @param method Method for aggregating across categories:
#'   {"none", "max", "sum"}. Default: "none".
#' @param weights Annotation category weights.
#' @return (n x 3) Numeric matrix without weighting, (n x 1) numeric matrix
#' with weighting.
#' @export
Aggregator <- function(
    anno,
    geno,
    drop_empty = TRUE,
    indicator = FALSE,
    method = "none",
    weights = DEFAULT_WEIGHTS
) {

  # Sum to categories.
  bmv <- apply(geno[, anno == 0, drop = FALSE], 1, sum)
  dmv <- apply(geno[, anno == 1, drop = FALSE], 1, sum)
  ptv <- apply(geno[, anno == 2, drop = FALSE], 1, sum)

  # Convert to indicators.
  if (indicator) {
    bmv <- 1 * (bmv > 0)
    dmv <- 1 * (dmv > 0)
    ptv <- 1 * (ptv > 0)
  }

  # If weights are null, output 3 column matrix.
  out <- cbind(bmv = bmv, dmv = dmv, ptv = ptv)

  # Apply weighting.
  out <- out %*% diag(weights)

  # Drop empty columns.
  if (drop_empty) {
    col_sums <- apply(out, 2, sum)
    out <- out[, col_sums > 0, drop = FALSE]

    # Case of no non-zero variant classes.
    if (ncol(out) == 0) {
      warning("No non-zero variant classes after weighting.")
      return(NA)
    }

  }

  if (method == "sum") {

    # Sum aggregation.
    out <- apply(out, 1, sum)

  } else if (method == "max") {

    # Max aggregation.
    out <- apply(out, 1, max)

  } else {

    # No aggregation (default).
    out <- out
  }

  return(as.matrix(out))
}


# ------------------------------------------------------------------------------


#' Allelic Series Burden Test
#'
#' Burden test with allelic series weights.
#'
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param geno (n x snps) genotype matrix.
#' @param pheno (n x 1) phenotype vector.
#' @param apply_int Apply rank-based inverse normal transform to the phenotype?
#'   Default: TRUE. Ignored if phenotype is binary.
#' @param covar (n x p) covariate matrix. Defaults to an (n x 1) intercept.
#' @param indicator Convert raw counts to indicators?
#' @param is_pheno_binary Is the phenotype binary? Default: FALSE.
#' @param method Method for aggregating across categories: {"none", "max",
#'   "sum"}. Default: "none".
#' @param score_test Run a score test? If FALSE, performs a Wald test.
#' @param weights (3 x 1) annotation category weights.
#' @return Numeric p-value.
#' @examples 
#' # Generate data.
#' data <- DGP(n = 1e3, snps = 1e2)
#' 
#' # Run the Allelic Series Burden Test.
#' # Note: the output is a scalar p-value.
#' results <- ASBT(
#'   anno = data$anno,
#'   geno = data$geno,
#'   pheno = data$pheno,
#'   covar = data$covar
#' )
#' @export
ASBT <- function(
  anno,
  geno,
  pheno,
  apply_int = TRUE,
  covar = NULL,
  indicator = FALSE,
  is_pheno_binary = FALSE,
  method = "none",
  score_test = FALSE,
  weights = DEFAULT_WEIGHTS
) {

  # Covariates.
  if (is.null(covar)) {
    covar <- rep(1, length(pheno))
  }

  # Data types.
  anno <- as.vector(anno)
  covar <- as.matrix(covar)
  geno <- as.matrix(geno)
  pheno <- as.vector(pheno)

  # Input check.
  CheckInputs(
    anno = anno,
    covar = covar,
    geno = geno,
    is_pheno_binary = is_pheno_binary,
    pheno = pheno,
    weights = weights
  )

  # Aggregate genotypes.
  agg_geno <- Aggregator(
    anno = anno,
    geno = geno,
    drop_empty = TRUE,
    indicator = indicator,
    method = method,
    weights = weights
  )

  # Case of no non-zero variant classes.
  if (!is.numeric(agg_geno)) {
    return(NA)
  }

  # Rank-normal phenotype.
  if (!is_pheno_binary & apply_int) {
    y <- RNOmni::RankNorm(pheno)
  } else {
    y <- pheno
  }

  # Association test.
  if (is_pheno_binary) {
    pval <- LogisticLRT(y = y, g = agg_geno, x = covar)
  } else {
    pval <- LinearTest(y = y, g = agg_geno, x = covar, score_test = score_test)
  }
  
  return(pval)
}


# ------------------------------------------------------------------------------

#' Allelic Series SKAT Test
#'
#' Sequence kernel association test (SKAT) with allelic series weights.
#'
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param geno (n x snps) genotype matrix.
#' @param pheno (n x 1) phenotype vector.
#' @param apply_int Apply rank-based inverse normal transform to the phenotype?
#'   Default: TRUE. Ignored if phenotype is binary.
#' @param covar (n x p) covariate matrix. Defaults to an (n x 1) intercept.
#' @param is_pheno_binary Is the phenotype binary? Default: FALSE.
#' @param return_null_model Return the null model in addition to the p-value?
#'   Useful if running additional SKAT tests. Default: FALSE.
#' @param weights (3 x 1) annotation category weights.
#' @return If `return_null_model`, a list containing the p-value and the
#'   SKAT null model. Otherwise, a numeric p-value.
#' @examples
#' # Generate data.
#' data <- DGP(n = 1e3, snps = 1e2)
#' 
#' # Run the Allelic Series SKAT Test.
#' # Note: the output is a scalar p-value.
#' results <- ASKAT(
#'   anno = data$anno,
#'   geno = data$geno,
#'   pheno = data$pheno,
#'   covar = data$covar
#' )
#' @export
ASKAT <- function(
  anno,
  geno,
  pheno,
  apply_int = TRUE,
  covar = NULL,
  is_pheno_binary = FALSE,
  return_null_model = FALSE,
  weights = DEFAULT_WEIGHTS
) {

  # Covariates.
  if (is.null(covar)) {
    covar <- rep(1, length(pheno))
  }

  # Data types.
  anno <- as.vector(anno)
  covar <- as.matrix(covar)
  geno <- as.matrix(geno)
  pheno <- as.vector(pheno)

  # Input check.
  CheckInputs(
    anno = anno,
    covar = covar,
    geno = geno,
    is_pheno_binary = is_pheno_binary,
    pheno = pheno,
    weights = weights
  )

  # Alternate allele frequencies.
  aaf <- apply(geno, 2, mean) / 2

  # Drop empty genotypes.
  is_empty <- (aaf == 0 | aaf == 1)

  anno <- anno[!is_empty]
  geno <- geno[, !is_empty, drop = FALSE]
  aaf <- aaf[!is_empty]

  # SKAT weight.
  w <- rep(0, length(anno))
  for (i in 0:2) {
    w[anno == i] <- weights[i + 1]
  }
  v <- aaf * (1 - aaf)
  skat_weights <- sqrt(w / v)

  # Case of no non-zero weights.
  if (all(skat_weights == 0)) {
    warning("No non-zero variant classes after weighting.")
    return(NA)
  }

  # Rank-normal phenotype.
  if (!is_pheno_binary & apply_int) {
    y <- RNOmni::RankNorm(pheno)
  } else {
    y <- pheno
  }

  # SKAT null model.
  if (is_pheno_binary) {
    null <- SKAT::SKAT_Null_Model(y ~ covar, out_type = "D")
  } else {
    null <- SKAT::SKAT_Null_Model(y ~ covar, out_type = "C")
  }

  # SKAT test.
  skato_test <- SKAT::SKAT(
    geno,
    null,
    method = "SKATO",
    weights = skat_weights
  )

  # Output.
  if (return_null_model) {
    out <- list(
      p = skato_test$p.value,
      null = null
    )
  } else {
    out <- skato_test$p.value
  }

  return(out)
}


# -----------------------------------------------------------------------------


#' COding-variant Allelic Series Test
#'
#' Main allelic series test. Performs both Burden and SKAT type tests, then
#' combines the results to calculate an omnibus p-value.
#'
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param geno (n x snps) genotype matrix.
#' @param pheno (n x 1) phenotype vector.
#' @param apply_int Apply rank-based inverse normal transform to the phenotype?
#'   Default: TRUE. Ignored if phenotype is binary.
#' @param covar (n x p) covariate matrix. Defaults to an (n x 1) intercept.
#' @param include_orig_skato_all Include the original version of SKAT-O applied
#' to all variants in the omnibus test? Default: FALSE.
#' @param include_orig_skato_ptv Include the original version of SKAT-O applied
#' to PTV variants only in the omnibus test? Default: FALSE.
#' @param is_pheno_binary Is the phenotype binary? Default: FALSE.
#' @param return_omni_only Return only the omnibus p-value? Default: FALSE.
#' @param score_test Use a score test for burden analysis? If FALSE, uses a 
#'   Wald test.
#' @param weights (3 x 1) annotation category weights.
#' @return Numeric p-value.
#' @examples 
#' # Generate data.
#' data <- DGP(n = 1e3, snps = 1e2)
#' 
#' # Run the COding-variant Allelic Series Test.
#' results <- COAST(
#'   anno = data$anno,
#'   geno = data$geno,
#'   pheno = data$pheno,
#'   covar = data$covar
#' )
#' show(results)
#' 
#' @export
COAST <- function(
  anno,
  geno,
  pheno,
  apply_int = TRUE,
  covar = NULL,
  include_orig_skato_all = FALSE,
  include_orig_skato_ptv = FALSE,
  is_pheno_binary = FALSE,
  return_omni_only = FALSE,
  score_test = FALSE,
  weights = DEFAULT_WEIGHTS
) {

  # Covariates.
  if (is.null(covar)) {
    covar <- rep(1, length(pheno))
  }

  # Data types.
  anno <- as.vector(anno)
  covar <- as.matrix(covar)
  geno <- as.matrix(geno)
  pheno <- as.vector(pheno)

  # Input check.
  CheckInputs(
    anno = anno,
    covar = covar,
    geno = geno,
    is_pheno_binary = is_pheno_binary,
    pheno = pheno,
    weights = weights
  )

  BurdenWrap <- function(...) {
    out <- ASBT(
      anno = anno,
      covar = covar,
      geno = geno,
      pheno = pheno,
      apply_int = apply_int,
      is_pheno_binary = is_pheno_binary,
      score_test = score_test,
      ...
    )
    return(out)
  }

  # Burden tests.
  p_count <- BurdenWrap(
    indicator = FALSE, method = "none", weights = c(1, 1, 1))

  # Case of no non-zero variant classes.
  if (is.na(p_count)) {return(NA)}

  p_ind <- BurdenWrap(
    indicator = TRUE, method = "none", weights = c(1, 1, 1))

  p_max_count <- BurdenWrap(
    indicator = FALSE, method = "max", weights = weights)

  p_max_ind <- BurdenWrap(
    indicator = TRUE, method = "max", weights = weights)

  p_sum_count <- BurdenWrap(
    indicator = FALSE, method = "sum", weights = weights)

  p_sum_ind <- BurdenWrap(
    indicator = TRUE, method = "sum", weights = weights)

  # Collect p-values.
  p_burden <- c(
    p_count = p_count,
    p_ind = p_ind,
    p_max_count = p_max_count,
    p_max_ind = p_max_ind,
    p_sum_count = p_sum_count,
    p_sum_ind = p_sum_ind
  )

  # SKAT tests.
  skat_list <- ASKAT(
    anno = anno,
    covar = covar,
    geno = geno,
    pheno = pheno,
    is_pheno_binary = is_pheno_binary,
    return_null_model = TRUE,
    weights = weights
  )
  p_skat <- c(p_allelic_skat = skat_list$p)
  null <- skat_list$null

  # Standard SKAT-O all.
  if (include_orig_skato_all) {
    orig_skato_all <- SKAT::SKAT(geno, null, method = "SKATO")
    p_skat <- c(p_skat, p_orig_skat_all = orig_skato_all$p.value)
  }

  # Standard SKAT-O PTV.
  if (include_orig_skato_ptv) {
    ptv <- geno[, anno == 2, drop = FALSE]
    if (ncol(ptv) > 0) {
      orig_skato_ptv <- SKAT::SKAT(ptv, null, method = "SKATO")
      p_skat <- c(p_skat, p_orig_skat_ptv = orig_skato_ptv$p.value)
    }
  }

  # Omnibus.
  n_burden <- length(p_burden)
  n_skat <- length(p_skat)

  p_val = c(p_burden, p_skat)
  omni_weights <- c(rep(1, n_burden), rep(n_burden / n_skat, n_skat))

  p_omni <- RNOmni::OmniP(p = p_val, w = omni_weights)

  # Output.
  if (return_omni_only) {
    out <- c(p_omni = p_omni)
  } else {
    out <- c(p_val, p_omni = p_omni)
  }
  return(out)
}
