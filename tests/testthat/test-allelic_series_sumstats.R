test_that("Check COAST on sumstats.", {
  
  # Null setting.
  withr::local_seed(123)
  data <- DGP(n = 1e3, prop_causal = 0)
  sumstats <- CalcSumstats(data = data)

  results <- COASTSS(
    anno = sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se,
    ld = sumstats$ld,
    maf = sumstats$maf
  )
  
  expect_true(all(results@Pvals$pval >= 0.05))
  
  # Alternative setting.
  withr::local_seed(123)
  data <- DGP(n = 1e3, prop_causal = 0.5)
  sumstats <- CalcSumstats(data = data)
  
  results <- COASTSS(
    anno = sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se,
    ld = sumstats$ld,
    maf = sumstats$maf
  )
  
  expect_true(all(results@Pvals$pval < 0.05))
  
})


# ------------------------------------------------------------------------------


test_that("Check COAST with missing inputs.", {
  
  # Null setting.
  withr::local_seed(123)
  data <- DGP(n = 1e3, prop_causal = 0)
  sumstats <- CalcSumstats(data = data)
  
  results <- suppressWarnings({COASTSS(
    anno = sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se
  )})
  
  expect_true(all(results@Pvals$pval >= 0.05))
  
  # Alternative setting.
  withr::local_seed(123)
  data <- DGP(n = 1e3, prop_causal = 0.5)
  sumstats <- CalcSumstats(data = data)
  
  results <- suppressWarnings({COASTSS(
    anno = sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se
  )})
  
  expect_true(all(results@Pvals$pval < 0.05))
  
})


