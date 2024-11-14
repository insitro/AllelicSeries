test_that("Overall check for binary phenotype.", {
  
  # Note: binary phenotypes require a larger sample size.
  withr::local_seed(101)
  n <- 100
  anno <- rep(c(1, 2, 3), each = 5)
  geno <- replicate(15, stats::rbinom(n, size = 2, prob = 0.1))
  eta <- -1 + geno %*% rep(c(0, 1, 2), each = 5)
  null_pheno <- stats::rbinom(n, size = 1, prob = 0.5)
  alt_pheno <- stats::rbinom(n, size = 1, prob = 1 / (1 + exp(-eta)))
  
  # Note: the invisible(capture.output({})) wrapper is used because SKAT
  # prints a small sample size warning.
  WrapCOAST <- function(pheno) {
    invisible(capture.output(
      out <- COAST(
        anno = anno,
        geno = geno,
        pheno = pheno,
        is_pheno_binary = TRUE,
        include_orig_skato_all = TRUE,
        include_orig_skato_ptv = TRUE
      )
    ))
    return(out)
  }
  
  # Null phenotype.
  p_omni <- WrapCOAST(null_pheno)
  pvals <- p_omni@Pvals
  p_omni <- as.numeric(pvals$pval[pvals$test == "omni"])
  expect_gt(p_omni, 0.05)
  
  # Alternative phenotype.
  p_omni <- WrapCOAST(alt_pheno)
  pvals <- p_omni@Pvals
  p_omni <- as.numeric(pvals$pval[pvals$test == "omni"])
  expect_equal(p_omni, 0.0, tolerance = 0.005)
  
})
  