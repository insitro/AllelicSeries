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


# ------------------------------------------------------------------------------

test_that("Check cases where not all annotaiton categories are provided.", {
  
  withr::local_seed(102)
  data <- DGP(n = 1e3, prop_causal = 0)

  ExcludeAnno <- function(data, anno) {
    key <- (data$anno != anno)
    out <- list(
      anno = data$anno[key],
      covar = data$covar,
      geno = data$geno[, key],
      pheno = data$pheno,
      type = data$type
    )
    return(out)
  }
  
  RunCOAST <- function(sumstats) {
    result <- COASTSS(
      anno = sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      ld = sumstats$ld,
      maf = sumstats$maf
    )
    pvals <- result@Pvals
    return(pvals)
  }
  
  # Exclude category 0.
  data0 <- ExcludeAnno(data, anno = 0)
  sumstats0 <- CalcSumstats(data = data0)
  result0 <- RunCOAST(sumstats0)
  expect_true(all(result0$pval >= 0.05))
  
  # Exclude category 1.
  data1 <- ExcludeAnno(data, anno = 1)
  sumstats1 <- CalcSumstats(data = data1)
  result1 <- RunCOAST(sumstats1)
  expect_true(all(result1$pval >= 0.05))
  
  # Exclude category 2.
  data2 <- ExcludeAnno(data, anno = 2)
  sumstats2 <- CalcSumstats(data = data2)
  result2 <- RunCOAST(sumstats2)
  expect_true(all(result2$pval >= 0.05))
  
  # Exclude two categories.
  data12 <- ExcludeAnno(data1, anno = 2)
  expect_true(all(data12$anno == 0))
  sumstats12 <- CalcSumstats(data = data12)
  result12 <- RunCOAST(sumstats12)
  expect_true(all(result12$pval >= 0.05))
  
})