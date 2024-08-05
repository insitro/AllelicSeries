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


# ------------------------------------------------------------------------------

test_that("Check genomic inflation factor.", {
  
  # Null setting.
  withr::local_seed(123)
  data <- DGP(n = 1e3, prop_causal = 0)
  sumstats <- CalcSumstats(data = data)
  
  uncorrected <- COASTSS(
    anno = sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se,
    ld = sumstats$ld,
    maf = sumstats$maf
  )
  p_uncorrected <- uncorrected@Pvals
  
  corrected <- suppressWarnings({COASTSS(
      anno = sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      lambda = c(1.2, 0.8, 1.0),
      ld = sumstats$ld,
      maf = sumstats$maf
    )
  })
  p_corrected <- corrected@Pvals

  expect_gt(p_corrected$pval[1], p_uncorrected$pval[1])  
  expect_equal(p_corrected$pval[2], p_uncorrected$pval[2])
  expect_equal(p_corrected$pval[3], p_uncorrected$pval[3])
  expect_gt(p_corrected$pval[4], p_uncorrected$pval[4])  
   
})


# ------------------------------------------------------------------------------

test_that("Check case of rank-deficient LD.", {
  
  withr::local_seed(101)
  data <- DGP(n = 1e2, prop_causal = 0)
  sumstats <- CalcSumstats(data = data)
  
  # The LD matrix is singular.
  expect_false(isPD(sumstats$ld))
  
  # Without epsilon, matrix inversion will fail.
  expect_error(
    COASTSS(
      anno = sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      check = FALSE,
      eps = 0,
      ld = sumstats$ld,
      maf = sumstats$maf
    )
  )
  
  # With epsilon, the test runs.
  results <- suppressWarnings(
    COASTSS(
      anno = sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      eps = 1e-4,
      ld = sumstats$ld,
      maf = sumstats$maf
    )
  )
  expect_true(all(results@Pvals$pval > 0.05))
  
  
})

