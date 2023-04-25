test_that("Check residual variance calculation.", {
  
  y <- c(-1, 1, 0, 1)
  x <- as.matrix(c(1, 1, 1, 1))
  
  # Observed variance. 
  obs_var <- ResidVar(y = y, X = x)
  
  # Expected variance.
  fit <- lm(y ~ 0 + x)
  fit_summary <- summary(fit)
  exp_var <- (fit_summary$sigma)^2
  
  expect_equal(obs_var, exp_var, ignore_attr = TRUE)  
  
})


test_that("Check score test calculation.", {
  
  y <- c(-1, 1, 0, 1)
  x <- as.matrix(c(1, 1, 1, 1))
  g <- as.matrix(c(2, 1, 0, 2))
  v <- 1
  
  # Observed score. 
  obs_score <- Score(y = y, G = g, X = x, v = v)
  
  # Expected score. 
  ey <- resid(lm(y ~ x))
  eg <- resid(lm(g ~ x))
  s <- t(g) %*% ey
  gpg <- t(eg) %*% eg
  exp_score <- as.numeric(t(s) %*% solve(gpg, s) / v)

  expect_equal(obs_score, exp_score, ignore_attr = TRUE)  
  
})


test_that("Overall check of score test.", {
  
  withr::local_seed(1010)
  
  # Null phenotype.
  data <- DGP(prop_causal = 0.0)
  p_omni_null <- COAST(
    anno = data$anno,
    geno = data$geno,
    pheno = data$pheno,
    include_orig_skato_all = FALSE,
    include_orig_skato_ptv = FALSE,
    score_test = TRUE
  )
  expect_gt(p_omni_null["p_omni"], 0.05)
  
  
  # Associated phenotype.
  data <- DGP(prop_causal = 1.0)
  p_omni_alt <- COAST(
    anno = data$anno,
    geno = data$geno,
    pheno = data$pheno,
    include_orig_skato_all = FALSE,
    include_orig_skato_ptv = FALSE,
    score_test = TRUE
  )
  expect_lt(p_omni_alt["p_omni"], 0.05)
  
})