# Note: Code modified from <https://github.com/leelabsg/SKAT>.
# Updated: 2024-08-02

#' Calculate P-value via Liu's Method
#' @param eigenvals Vector of eigenvalues.
#' @param q SKAT Q statistic.
#' @return Data.frame of test results.
#' @noRd
LiuMethod <- function(eigenvals, q) {
  
  # Moments.
  m1 <- sum(eigenvals)
  m2 <- sum(eigenvals^2)
  m3 <- sum(eigenvals^3)
  m4 <- sum(eigenvals^4)
  
  s1 <- m3 / m2^(3 / 2)
  s2 <- m4 / m2^2
  
  # Calculate degrees of freedom.
  if (s1^2 > s2) {
    a <- 1 / (s1 - sqrt(s1^2 - s2))
    d <- s1 * a^3 - a^2
    df <- a^2 - 2 * d
  } else {
    df <- 1 / s2
  }
  
  # Calculate p-value.
  mu_q <- m1
  sigma_q <- sqrt(2 * m2)
  q_norm <- (q - mu_q) / sigma_q * sqrt(2 * df) + df
  
  # Output.
  out <- data.frame(
    qstat = q,
    df = df,
    mu_q = mu_q,
    sigma_q = sigma_q,
    p_liu = stats::pchisq(q_norm, df = df, lower.tail = FALSE)
  )
  return(out)
}


#' Davies P-value
#' @param lambda Vector of eigenvalues.
#' @param q SKAT Q statistic.
#' @param acc Accuracy.
#' @param p_default Fallback p-value.
#' @return Numeric p-value.
#' @noRd
DaviesP <- function(
    lambda, q, p_default = NA, acc = 1e-8) {
  
  # Davies' method.
  davies_output <- suppressWarnings({
    tryCatch(
      CompQuadForm::davies(q = q, lambda = lambda, acc = acc),
      error = function(cond) {return(NULL)}
    )
  })
  
  # Case of Davies' method failure.
  if (is.null(davies_output) || 
      davies_output$ifault != 0 || 
      davies_output$Qq <= 0) {
    return(p_default)
  }

  return(davies_output$Qq)
}


#' Calculate Results for each Rho
#' @param lambdas List of eigenvalues.
#' @param q_stats SKAT Q statistics.
#' @param rhos Rho values.
#' @noRd
PerRhoResults <- function(
    lambdas,
    qstats,
    rhos
) {
  
  # Number of rhos considered.
  n_r <- length(rhos)
  
  results <- lapply(seq_len(n_r), function(i) {
    
    # Liu's method.
    out <- LiuMethod(lambdas[[i]], qstats[i])
    
    # Davies' method.
    out$p_davies <- DaviesP(
      lambda = lambdas[[i]],
      q = qstats[[i]],
      p_default = out$p_liu
    )
    
    return(out)
  })
  results <- do.call(rbind, results)
  
  # Update Q-statistics.
  pmin <- min(results$p_davies)
  results$qstat <- sapply(seq_len(n_r), function(i) {
    mu_q <- results$mu_q[i]
    sigma_q <- results$sigma_q[i]
    df <- results$df[i]
    q <- stats::qchisq(pmin, df = df, lower.tail = FALSE)
    q <- (q - df) / sqrt(2 * df) * sigma_q + mu_q
  })
  
  # Output.
  results$rho <- rhos
  return(results)
}


# ------------------------------------------------------------------------------

#' Davies Integrand
#' @param x Scalar variable of integration.
#' @param opt_params SKAT optimal parameters.
#' @param results Output of \code{\link{PerRhoResults}}.
#' @return Numeric p-value
#' @noRd 
DaviesIntegrand <- function(x, opt_params, results){
  
  x_tau <- as.numeric(x * opt_params$tau)
  qstats <- (results$qstat - x_tau) / (1 - results$rho) 
  min_q <- min(qstats)
  
  if (min_q > sum(opt_params$lambda) * 1e4) {
    alpha <- 0
  } else {
    mu_q <- opt_params$mu_q
    var_q <- opt_params$var_q
    var_e <- opt_params$var_remain
    se <- sqrt(var_q - var_e) / sqrt(var_q)
    q_stat <- (min_q - mu_q) * se + mu_q
    alpha <- DaviesP(lambda = opt_params$lambda, q = q_stat)
  }
  alpha <- min(alpha, 1)
  out <- (1 - alpha) * stats::dchisq(x, df = 1)
  return(out)
}


#' Liu Integrand
#' @param x Scalar variable of integration.
#' @param opt_params SKAT optimal parameters.
#' @param results Output of \code{\link{PerRhoResults}}.
#' @return Numeric p-value
#' @noRd
LiuIntegrand <- function(x, opt_params, results){
 
  x_tau <- as.numeric(x * opt_params$tau)
  qstats <- (results$qstat - x_tau) / (1 - results$rho) 
  min_q <- min(qstats)
  mu_q <- opt_params$mu_q
  var_q <- opt_params$var_q
  df <- opt_params$df
  
  q_stat <- (min_q - mu_q) / sqrt(var_q) * sqrt(2 * df) + df
  out <- stats::pchisq(q_stat, df = df) * stats::dchisq(x, df = 1) 
  return(out)
}


# ------------------------------------------------------------------------------

#' Optimal P-value
#' @param opt_params SKAT optimal parameters.
#' @param results Output of \code{\link{PerRhoResults}}.
#' @return Numeric p-value.
#' @noRd 
OptimalPval <-function(opt_params, results){
  
  # Try Davies' method first.
  aux <- function(x) {DaviesIntegrand(x, opt_params, results)}
  aux <- Vectorize(aux)
  res <- tryCatch(
    stats::integrate(aux, lower = 0, upper = 40, subdivisions = 2e3),
    error = function(cond) {return(NULL)}
  )
  
  if (!is.null(res) && res$value < 1) {
    pval <- 1 - res$value
  } else {
    
    # Fallback to Liu's method.
    aux <- function(x) {LiuIntegrand(x, opt_params, results)}
    aux <- Vectorize(aux)
    res <- stats::integrate(aux, lower = 0, upper = 40, subdivisions = 2e3)
    pval <- 1 - res$value
  }

  pmin <- min(results$p_davies)
  n_r <- nrow(results)
  if (pmin * n_r < pval) {
    pval <- pmin * n_r
  }

  return(pval)
}

