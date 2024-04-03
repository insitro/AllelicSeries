# Purpose: Utility functions.
# Updated: 2024-04-03

#' Count Alleles
#'
#' Count the number of non-zero alleles bearing each variant annotation.
#' 
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param geno (n x snps) genotype matrix.
#' @param count_carriers If true, counts the number of carriers rather than
#'   the number of alleles.
#' @param min_mac Minimum minor allele count for inclusion. Default: 0. 
#' @return 3 x 1 numeric vector with the counts of BMVs, DMVs, and PTVs.
#' @export
CountAlleles <- function(
    anno, 
    geno, 
    count_carriers = FALSE,
    min_mac = 0
) {
  
  # Minor allele count filtering.
  if (min_mac > 0) {
    mac <- apply(geno, 2, sum)
    keep <- (mac >= min_mac)
    
    anno <- anno[keep]
    geno <- geno[, keep, drop = FALSE]
  }
  
  # Category counts.
  if (count_carriers) {
    
    # Count carriers. 
    n_bmv <- sum(apply(geno[, anno == 0, drop = FALSE], 1, sum) > 0)
    n_dmv <- sum(apply(geno[, anno == 1, drop = FALSE], 1, sum) > 0)
    n_ptv <- sum(apply(geno[, anno == 2, drop = FALSE], 1, sum) > 0)
    
  } else {
    
    # Count alleles.
    n_bmv <- sum(geno[, anno == 0])
    n_dmv <- sum(geno[, anno == 1])
    n_ptv <- sum(geno[, anno == 2])
    
  }
  
  # Output.
  out <- c(n_bmv = n_bmv, n_dmv = n_dmv, n_ptv = n_ptv)
  return(out)
}


#' Linear Association Test
#' 
#' @param y (n x 1) continuous phenotype vector
#' @param g (n x p) genotype matrix.
#' @param x (n x q) covariate matrix.
#' @param score_test Run a score test? If FALSE, performs a Wald test. 
#' @return Numeric p-value.
#' @noRd
LinearTest <- function(y, g, x, score_test = FALSE) {
  
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
  return(pval)
  
}


#' Likelihood Ratio Test
#' 
#' @param y (n x 1) binary (0/1) phenotype vector
#' @param g (n x p) genotype matrix.
#' @param x (n x q) covariate matrix.
#' @return Numeric p-value.
#' @noRd
LogisticLRT <- function(y, g, x) {
  
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
  return(pval)
  
}
