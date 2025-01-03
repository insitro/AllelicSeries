# Purpose: Utility functions.
# Updated: 2025-01-02

#' Collapse Variants
#' 
#' Collapse variants with minor allele counts below the `min_mac` threshold
#' into an aggregated variant, separately within each variant category. Note 
#' that the ordering of the variants will change, and that collapsing does not
#' guarantee that the resulting aggregate aggregate variant will itself have a 
#' MAC greater than or equal to `min_mac`.
#' 
#' @param anno (snps x 1) annotation vector with integer values in 1 through
#'   the number of annotation categories L.
#' @param geno (n x snps) genotype matrix.
#' @param min_mac Minimum minor allele count (MAC) for retention as an 
#'   individual variant. Variants with a MAC strictly less than the minimum
#'   MAC will be collapsed into an aggregated variant, separately within
#'   each annotation category.
#' @return List containing the collapsed genotypes `geno` and corresponding 
#'   annotations `anno`, plus a data.frame `vars` specifying which variants
#'   were collapsed within each annotation category.
#' @export
CollapseGeno <- function(
    anno,
    geno,
    min_mac = 11
) {
  
  # Check for missing values.
  if (any(is.na(anno), is.na(geno))) {
    stop("Missing values should be addressed before aggregation.")
  }
  
  # Filter out variants with a MAC of zero.
  macs <- apply(geno, 2, sum)
  if (any(macs == 0)) {
    warning("Variants with a MAC of zero are present. These will be dropped.")
    key <- (macs > 0) 
    anno <- anno[key]
    geno <- geno[, key, drop = FALSE]
    macs <- macs[key]
  }
  if (ncol(geno) == 0) {return(NULL)}
  
  # Loop over annotations.
  annos_out <- c()
  collapsed <- array(dim = c(nrow(geno), 0))
  unique_annos <- sort(unique(anno))
  vars_agg <- data.frame(anno = unique_annos, vars = NA)
  
  for(l in unique_annos) {
    
    # Subset to current annotation.
    geno_l <- geno[, anno == l, drop = FALSE]
    macs_l <- macs[anno == l]
    
    # Partition into individual variants and the aggregated variant.
    keep <- (macs_l >= min_mac)
    indiv <- geno_l[, keep, drop = FALSE]
    to_agg <- geno_l[, !keep, drop = FALSE]
    if (ncol(to_agg) > 0) {
      agg <- apply(to_agg, 1, sum)
    } else {
      agg <- NULL
    }
    
    # Variants aggregated.
    vars_agg$vars[vars_agg$anno == l] <- paste(colnames(to_agg), collapse = ";")
    
    # Collapsed genotypes.
    collapsed_l <- cbind(indiv, agg)
    if (!is.null(agg)) {
      colnames(collapsed_l)[ncol(indiv) + 1] <- glue::glue("agg{l}")
    }
    collapsed <- cbind(collapsed, collapsed_l)
    
    # Annotations.
    annos_out <- c(annos_out, rep(l, times = ncol(collapsed_l)))
  }
  
  # Output.
  out <- list(
    anno = annos_out,
    geno = collapsed,
    vars = vars_agg
  )
  return(out)
}


#' Genomic Control
#' 
#' @param lambda Genomic inflation factor.
#' @param pval Numeric p-value.
#' @param df Degrees of freedom. Should not require modification in most cases.
#' @return Corrected p-value.
#' @export 
GenomicControl <- function(lambda, pval, df = 1) {
  if (lambda == 1 || is.na(pval)) {
    return(pval)
  }
  chi2 <- stats::qchisq(p = pval, df = df, lower.tail = FALSE) 
  chi2_corrected <- chi2 / lambda
  pval_corrected <- stats::pchisq(
    q = chi2_corrected, df = df, lower.tail = FALSE)
  return(pval_corrected)
}


#' Linear Association Test
#' 
#' @param y (n x 1) continuous phenotype vector
#' @param g (n x p) genotype matrix.
#' @param x (n x q) covariate matrix.
#' @param return_beta Return the estimated effect size? Default: FALSE.
#' @param score_test Run a score test? If FALSE, performs a Wald test. 
#' @return If `return_beta = TRUE`, a list of including the effect size
#'   data.frame "betas" and the p-value "pval". If `return_beta = FALSE`, 
#'   a numeric p-value.
#' @noRd
LinearTest <- function(y, g, x, return_beta = FALSE, score_test = FALSE) {
  
  # Estimate residual variance.
  if (score_test) {
    resid_var <- ResidVar(y = y, X = x)
  } else {
    resid_var <- ResidVar(y = y, X = cbind(g, x))
  }
  
  # Calculate test statistic.
  stat <- Score(y = y, G = g, X = x, v = resid_var)
  df <- ncol(g)
  pval <- stats::pchisq(q = stat, df = df, lower.tail = FALSE)
  
  # Return p-value only if effect sizes not requested.
  if (!return_beta) {return(pval)}
  
  # Calculate effect sizes.
  fit <- OLS(y = y, X = cbind(g, x))
  n_beta <- ncol(g)
  betas <- data.frame(
    beta = fit$beta[1:n_beta],
    se = fit$se[1:n_beta]
  )
  
  # Return effect sizes and p-values.
  out <- list(
    betas = betas,
    pval = pval
  )
  return(out)
}


#' Likelihood Ratio Test
#' 
#' @param y (n x 1) binary (0/1) phenotype vector
#' @param g (n x p) genotype matrix.
#' @param x (n x q) covariate matrix.
#' @param return_beta Return the estimated effect size? Default: FALSE.
#' @return If `return_beta = TRUE`, a list of including the effect size
#'   data.frame "betas" and the p-value "pval". If `return_beta = FALSE`, 
#'   a numeric p-value.
#' @noRd
LogisticLRT <- function(y, g, x, return_beta = FALSE) {
  
  # Full model.
  fit_full <- stats::glm(
      y ~ 0 + g + x, family = stats::binomial(link = "logit"))
  
  # Null model.
  fit_null <- stats::glm(
    y ~ 0 + x, family = stats::binomial(link = "logit"))
  
  # LRT.
  lrt <- stats::anova(fit_null, fit_full)
  dev <- lrt$Deviance[2]
  df <- lrt$Df[2]
  pval <- stats::pchisq(q = dev, df = df, lower.tail = FALSE)
  
  # Return p-value only if effect sizes not requested.
  if (!return_beta) {return(pval)}
  
  # Calculate effect sizes.
  fit_summary <- summary(fit_full)
  n_beta <- ncol(g)
  betas <- data.frame(
    beta = fit_summary$coefficients[1:n_beta, 1],
    se = fit_summary$coefficients[1:n_beta, 2]
  )
  
  # Return effect sizes and p-values.
  out <- list(
    betas = betas,
    pval = pval
  )
  return(out)
}

