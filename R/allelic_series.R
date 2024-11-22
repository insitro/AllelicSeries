# Purpose: Allelic series test.
# Updated: 2024-11-20

#' Aggregator
#'
#' Aggregates genotypes within annotation categories. 
#' 
#' @section Notes:
#' * Ensure the length of the `weights` vector matches the total number of 
#'   annotation categories.
#' * The `weights` essentially scales the minor allele count in the `l`th
#'   category by `weights[l]`.
#'
#' @param anno (snps x 1) annotation vector with integer values in 1 through
#'   the number of annotation categories L.
#' @param geno (n x snps) genotype matrix.
#' @param drop_empty Drop empty columns? Default: TRUE.
#' @param indicator Convert raw counts to indicators? Default: FALSE.
#' @param method Method for aggregating across categories:
#'   ("none", "max", "sum"). Default: "none".
#' @param min_mac Minimum minor allele count for inclusion. Default: 0. 
#' @param weights (L x 1) vector of annotation category weights. Note that the
#'   number of annotation categories L is inferred from the length of `weights`.
#' @return (n x L) Numeric matrix without weighting, (n x 1) numeric matrix
#' with weighting.
#' @export
Aggregator <- function(
    anno,
    geno,
    drop_empty = TRUE,
    indicator = FALSE,
    method = "none",
    min_mac = 0,
    weights = c(1, 2, 3)
) {
  
  # Ensure annotations are one-indexed.
  anno <- RelevelAnno(anno)
  
  # Minor allele count filtering.
  if (min_mac > 0) {
    mac <- apply(geno, 2, sum)
    keep <- (mac >= min_mac)
    
    anno <- anno[keep]
    geno <- geno[, keep, drop = FALSE]
  }

  # Sum to categories.
  n_anno <- length(weights)
  out <- lapply(seq_len(n_anno), function(l) {
    agg <- apply(geno[, anno == l, drop = FALSE], 1, sum)
    if (indicator) {agg <- 1 * (agg > 0)}
    return(agg)
  })
  out <- do.call(cbind, out)
  colnames(out) <- paste0("a", seq_len(n_anno))
  
  stopifnot(ncol(out) == n_anno)
  stopifnot(nrow(out) == nrow(geno))

  # Apply annotation category weights
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
#' @param anno (snps x 1) annotation vector with integer values in 1 through
#'   the number of annotation categories L.
#' @param geno (n x snps) genotype matrix.
#' @param pheno (n x 1) phenotype vector.
#' @param apply_int Apply rank-based inverse normal transform to the phenotype?
#'   Default: TRUE. Ignored if phenotype is binary.
#' @param covar (n x p) covariate matrix. Defaults to an (n x 1) intercept.
#' @param indicator Convert raw counts to indicators?
#' @param is_pheno_binary Is the phenotype binary? Default: FALSE.
#' @param method Method for aggregating across categories: ("none", "max",
#'   "sum"). Default: "none".
#' @param min_mac Minimum minor allele count for inclusion. Default: 0. 
#' @param return_beta Return the estimated effect size? Default: FALSE.
#' @param score_test Run a score test? If FALSE, performs a Wald test.
#' @param weights (L x 1) vector of annotation category weights. Note that the
#'   number of annotation categories L is inferred from the length of `weights`.
#' @return If `return_beta = TRUE`, a list of including the effect size
#'   data.frame "betas" and the p-value "pval". If `return_beta = FALSE`, 
#'   a numeric p-value.
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
  min_mac = 0,
  return_beta = FALSE,
  score_test = FALSE,
  weights = c(1, 2, 3)
) {
  
  # Ensure annotations are one-indexed.
  anno <- RelevelAnno(anno)

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
    min_mac = min_mac,
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
    out <- LogisticLRT(
      y = y, 
      g = agg_geno, 
      x = covar,
      return_beta = return_beta
    )
  } else {
    out <- LinearTest(
      y = y, 
      g = agg_geno, 
      x = covar, 
      return_beta = return_beta,
      score_test = score_test
    )
  }
  
  return(out)
}


# ------------------------------------------------------------------------------

#' Allelic Series SKAT Test
#'
#' Sequence kernel association test (SKAT) with allelic series weights.
#'
#' @param anno (snps x 1) annotation vector with integer values in 1 through
#'   the number of annotation categories L.
#' @param geno (n x snps) genotype matrix.
#' @param pheno (n x 1) phenotype vector.
#' @param apply_int Apply rank-based inverse normal transform to the phenotype?
#'   Default: TRUE. Ignored if phenotype is binary.
#' @param covar (n x p) covariate matrix. Defaults to an (n x 1) intercept.
#' @param is_pheno_binary Is the phenotype binary? Default: FALSE.
#' @param min_mac Minimum minor allele count for inclusion. Default: 0. 
#' @param return_null_model Return the null model in addition to the p-value?
#'   Useful if running additional SKAT tests. Default: FALSE.
#' @param weights (L x 1) vector of annotation category weights. Note that the
#'   number of annotation categories L is inferred from the length of `weights`.
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
  min_mac = 0,
  return_null_model = FALSE,
  weights = c(1, 2, 3)
) {
  
  # Ensure annotations are one-indexed.
  anno <- RelevelAnno(anno)

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

  # Minor allele count filtering.
  if (min_mac >= 0) {
    mac <- apply(geno, 2, sum)
    keep <- (mac > min_mac)
    
    anno <- anno[keep]
    geno <- geno[, keep, drop = FALSE]
  }
  
  # Alternate allele frequencies.
  aaf <- apply(geno, 2, mean) / 2

  # Drop empty genotypes.
  is_empty <- (aaf == 0 | aaf == 1)

  anno <- anno[!is_empty]
  geno <- geno[, !is_empty, drop = FALSE]
  aaf <- aaf[!is_empty]

  # SKAT weight.
  w <- rep(0, length(anno))
  n_anno <- length(weights)
  for (i in seq_len(n_anno)) {
    w[anno == i] <- weights[i]
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
#' @param anno (snps x 1) annotation vector with integer values in 1 through
#'   the number of annotation categories L.
#' @param geno (n x snps) genotype matrix.
#' @param pheno (n x 1) phenotype vector.
#' @param apply_int Apply rank-based inverse normal transform to the phenotype?
#'   Default: TRUE. Ignored if phenotype is binary.
#' @param covar (n x p) covariate matrix. Defaults to an (n x 1) intercept.
#' @param include_orig_skato_all Include the original version of SKAT-O applied
#'   to all variants in the omnibus test? Default: FALSE.
#' @param include_orig_skato_ptv Include the original version of SKAT-O applied
#'   to PTV variants only in the omnibus test? Default: FALSE.
#' @param is_pheno_binary Is the phenotype binary? Default: FALSE.
#' @param min_mac Minimum minor allele count for inclusion. Default: 0. 
#' @param ptv_anno Annotation of the PTV category, only required if 
#'   include_orig_skato_ptv is set to TRUE.
#' @param pval_weights Optional vector of relative weights for combining the
#'   component tests to perform the omnibus test. By default, 50% of weight is
#'   given to the 6 burden tests, and 50% to the 1 SKAT test. If specified, the
#'   weight vector should have length 7, and the length should be increased if
#'   either `include_orig_skato_all` or `include_orig_skato_ptv` is active.
#' @param return_omni_only Return only the omnibus p-value? Default: FALSE.
#' @param score_test Use a score test for burden analysis? If FALSE, uses a 
#'   Wald test.
#' @param weights (L x 1) vector of annotation category weights. Note that the
#'   number of annotation categories L is inferred from the length of `weights`.
#' @return An object of class `COAST` with slots for effect sizes, variant 
#'   counts, and p-values.
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
  min_mac = 0,
  ptv_anno = 3,
  pval_weights = NULL,
  return_omni_only = FALSE,
  score_test = FALSE,
  weights = c(1, 2, 3)
) {

  # Ensure annotations are one-indexed.
  anno <- RelevelAnno(anno)
  
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
      min_mac = min_mac,
      is_pheno_binary = is_pheno_binary,
      return_beta = TRUE,
      score_test = score_test,
      ...
    )
    return(out)
  }

  # Burden tests.
  n_anno <- length(weights)
  unif_weights <- rep(1, n_anno)
  
  burden_results <- list()
  burden_results$base <- BurdenWrap(
    indicator = FALSE, method = "none", weights = unif_weights)

  # Case of no non-zero variant classes.
  if (is.na(burden_results$base$pval)) {return(NA)}

  burden_results$ind <- BurdenWrap(
    indicator = TRUE, method = "none", weights = unif_weights)

  burden_results$max_count <- BurdenWrap(
    indicator = FALSE, method = "max", weights = weights)

  burden_results$max_ind <- BurdenWrap(
    indicator = TRUE, method = "max", weights = weights)

  burden_results$sum_count <- BurdenWrap(
    indicator = FALSE, method = "sum", weights = weights)

  burden_results$sum_ind <- BurdenWrap(
    indicator = TRUE, method = "sum", weights = weights)

  # Collect p-values.
  p_burden <- c(
    baseline = burden_results$base$pval,
    ind = burden_results$ind$pval,
    max_count = burden_results$max_count$pval,
    max_ind = burden_results$max_ind$pval,
    sum_count = burden_results$sum_count$pval,
    sum_ind = burden_results$sum_ind$pval
  )

  # SKAT tests.
  skat_list <- ASKAT(
    anno = anno,
    covar = covar,
    geno = geno,
    pheno = pheno,
    is_pheno_binary = is_pheno_binary,
    min_mac = min_mac,
    return_null_model = TRUE,
    weights = weights
  )
  p_skat <- c(allelic_skat = skat_list$p)
  null <- skat_list$null

  # Standard SKAT-O all.
  if (include_orig_skato_all) {
    orig_skato_all <- SKAT::SKAT(geno, null, method = "SKATO")
    p_skat <- c(p_skat, orig_skat_all = orig_skato_all$p.value)
  }

  # Standard SKAT-O PTV.
  if (include_orig_skato_ptv) {
    ptv <- geno[, anno == ptv_anno, drop = FALSE]
    if (ncol(ptv) > 0) {
      orig_skato_ptv <- SKAT::SKAT(ptv, null, method = "SKATO")
      p_skat <- c(p_skat, orig_skat_ptv = orig_skato_ptv$p.value)
    }
  }

  # Omnibus.
  n_burden <- length(p_burden)
  n_skat <- length(p_skat)

  pvals = c(p_burden, p_skat)
  omni_weights <- c(rep(1, n_burden), rep(n_burden / n_skat, n_skat))
  if (!is.null(pval_weights)) {
    len_pval_weights <- length(pval_weights)
    len_omni_weights <- length(omni_weights)
    if (len_pval_weights == len_omni_weights) {
      omni_weights <- pval_weights
    } else {
      msg <- paste0(
        glue::glue("pval_weights has length {len_pval_weights}, but length {len_omni_weights} is needed. "),
        "Default weights will be used instead."
      )
      warning(msg)
    }
  }

  p_omni <- RNOmni::OmniP(p = pvals, w = omni_weights)

  # Only omnibus p-value requested.
  if (return_omni_only) {
    out <- c(p_omni = p_omni)
    return(out)
  } 
  
  # Format p-values.
  pvals <- c(pvals, omni = p_omni)
  df_pvals <- data.frame(
    test = names(pvals),
    type = c(rep("burden", n_burden), rep("skat", n_skat), "omni"),
    pval = as.numeric(pvals)
  )
  
  # Variant counts.
  counts <- Counts(
    anno = anno, 
    geno = geno, 
    n_anno = n_anno,
    min_mac = min_mac
  )
  
  # Format betas.
  burden_labs <- names(burden_results)
  betas <- lapply(burden_labs, function(lab) {
    out <- burden_results[[lab]]$betas
    out$test <- lab
    out <- out[, c("test", "beta", "se")]
    return(out)
  })
  betas <- do.call(rbind, betas)
  
  # Output.
  out <- methods::new(
    Class = "COAST",
    Betas = betas,
    Counts = counts,
    Pvals = df_pvals
  )
  return(out)
}
