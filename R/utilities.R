# Purpose: Utility functions.
# Updated: 2022-09-22

#' Fit Logistic Regression
#'
#' Potential extension: write an estimation routine in C++ (to increase speed).
#'
#' @param y (n x 1) binary (0/1) phenotype vector.
#' @param x (n x q) design matrix including both the aggregated genotypes and
#'   the covariates.
#' @return List.
#' @noRd
GLM <- function(y, x) {
  fit <- stats::glm(y ~ 0 + x, family = stats::binomial(link = "logit"))
  fit_summary <- stats::summary.glm(fit)
  fit_coef <- fit_summary$coefficients
  out <- list(z = fit_coef[, 3])
  return(out)
}
