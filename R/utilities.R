# Purpose: Utility functions.
# Updated: 2024-04-03

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
