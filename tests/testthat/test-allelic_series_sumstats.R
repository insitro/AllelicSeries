test_that("Check COAST on sumstats.", {
  
  # Null setting.
  withr::local_seed(123)
  data <- DGP(n = 1e3, prop_causal = 0)
  sumstats <- CalcSumstats(data = data)

  results <- COASTSS(
    anno = sumstats$sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se,
    ld = sumstats$ld,
    maf = sumstats$sumstats$maf
  )
  
  expect_true(all(results@Pvals$pval >= 0.05))
  
  # Alternative setting.
  withr::local_seed(123)
  data <- DGP(n = 1e3, prop_causal = 0.5)
  sumstats <- CalcSumstats(data = data)
  
  results <- COASTSS(
    anno = sumstats$sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se,
    ld = sumstats$ld,
    maf = sumstats$sumstats$maf
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
    anno = sumstats$sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se
  )})
  
  expect_true(all(results@Pvals$pval >= 0.05))
  
  # Alternative setting.
  withr::local_seed(123)
  data <- DGP(n = 1e3, prop_causal = 0.5)
  sumstats <- CalcSumstats(data = data)
  
  results <- suppressWarnings({COASTSS(
    anno = sumstats$sumstats$anno,
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
    anno = sumstats$sumstats$anno,
    beta = sumstats$sumstats$beta,
    se = sumstats$sumstats$se,
    ld = sumstats$ld,
    maf = sumstats$sumstats$maf
  )
  p_uncorrected <- uncorrected@Pvals
  
  corrected <- suppressWarnings({COASTSS(
      anno = sumstats$sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      lambda = c(1.2, 0.8, 1.0),
      ld = sumstats$ld,
      maf = sumstats$sumstats$maf
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
  
  withr::local_seed(102)
  data <- DGP(n = 1e2, prop_causal = 0)
  sumstats <- CalcSumstats(data = data)
  ld <- sumstats$ld
  
  # Ensure the matrix is singular.
  ld[1, 0] <- ld[0, 1] <- 0
  
  # The LD matrix is singular.
  expect_false(isPD(sumstats$ld))
  
  # Test COAST runs even *without* epsilon.
  expect_error(
    COASTSS(
      anno = sumstats$sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      check = FALSE,
      eps = 0,
      ld = sumstats$ld,
      maf = sumstats$sumstats$maf
    ), 
    NA
  )
  
  # With epsilon, the test runs.
  results <- suppressWarnings(
    COASTSS(
      anno = sumstats$sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      eps = 1,
      ld = sumstats$ld,
      maf = sumstats$sumstats$maf
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
    result <- suppressWarnings(
      COASTSS(
        anno = sumstats$sumstats$anno,
        beta = sumstats$sumstats$beta,
        se = sumstats$sumstats$se,
        ld = sumstats$ld,
        maf = sumstats$sumstats$maf
      )
    )
    pvals <- result@Pvals
    return(pvals)
  }
  
  # Exclude category 1.
  data0 <- ExcludeAnno(data, anno = 1)
  sumstats0 <- CalcSumstats(data = data0)
  result0 <- RunCOAST(sumstats0)
  expect_true(all(result0$pval >= 0.05))
  
  # Exclude category 2.
  data1 <- ExcludeAnno(data, anno = 2)
  sumstats1 <- CalcSumstats(data = data1)
  result1 <- RunCOAST(sumstats1)
  expect_true(all(result1$pval >= 0.05))
  
  # Exclude category 3.
  data2 <- ExcludeAnno(data, anno = 3)
  sumstats2 <- CalcSumstats(data = data2)
  result2 <- RunCOAST(sumstats2)
  expect_true(all(result2$pval >= 0.05))
  
  # Exclude two categories.
  data12 <- ExcludeAnno(data1, anno = 3)
  expect_true(all(data12$anno == 1))
  sumstats12 <- CalcSumstats(data = data12)
  result12 <- RunCOAST(sumstats12)
  expect_true(all(result12$pval >= 0.05))
  
})

# ------------------------------------------------------------------------------

test_that("Check COASTSS runs with zero-based annotations.", {
  
  withr::local_seed(103)
  
  WrapCOASTSS <- function(data, weights) {
    sumstats <- CalcSumstats(data = data)
    test <- COASTSS(
      anno = sumstats$sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      ld = sumstats$ld,
      maf = sumstats$sumstats$maf
    )
    pvals <- test@Pvals
    p_omni <- pvals$pval[pvals$test == "omni"]
    return(p_omni)
  }
  
  # Null. 
  data <- DGP(prop_causal = 0)
  data$anno <- data$anno - 1
  suppressWarnings(
    expect_warning(p <- WrapCOASTSS(data))
  )
  expect_gt(p, 0.05)
  
  # Alternative. 
  data <- DGP(prop_causal = 1)
  data$anno <- data$anno - 1
  suppressWarnings(
    expect_warning(p <- WrapCOASTSS(data))
  )
  expect_lt(p, 0.05)
  
})

# ------------------------------------------------------------------------------

test_that("Check COASTSS with different numbers of categories.", {
  
  WrapCOASTSS <- function(data, weights) {
    sumstats <- CalcSumstats(data = data)
    test <- suppressWarnings({COASTSS(
      anno = sumstats$sumstats$anno,
      beta = sumstats$sumstats$beta,
      se = sumstats$sumstats$se,
      ld = sumstats$ld,
      maf = sumstats$sumstats$maf,
      weights = weights
    )})
    pvals <- test@Pvals
    p_omni <- pvals$pval[pvals$test == "omni"]
    return(p_omni)
  }
  
  withr::local_seed(101)
  
  # 2 categories.
  ## Null.
  data <- DGP(
    beta = c(0, 0),
    prop_anno = c(1, 1),
    weights = c(1, 1)
  )
  p <- WrapCOASTSS(data, weights = c(1, 2))
  expect_gt(p, 0.05)
  
  ## Alternative.
  data <- DGP(
    beta = c(1, 2),
    prop_anno = c(1, 1),
    weights = c(1, 1)
  )
  p <- WrapCOASTSS(data, weights = c(1, 2))
  expect_lt(p, 0.05)
  
  # 4 categories.
  ## Null.
  data <- DGP(
    beta = c(0, 0, 0, 0),
    prop_anno = c(1, 1, 1, 1),
    weights = c(1, 1, 1, 1)
  )
  p <- WrapCOASTSS(data, weights = c(1, 2, 3, 4))
  expect_gt(p, 0.05)
  
  ## Alternative.
  data <- DGP(
    beta = c(1, 2, 3, 4),
    prop_anno = c(1, 1, 1, 1),
    weights = c(1, 1, 1, 1)
  )
  p <- WrapCOASTSS(data, weights = c(1, 2, 3, 4))
  expect_lt(p, 0.05)
  
})

# ------------------------------------------------------------------------------

# Slow.
test_that("Test effect size estimation from sumstats.", {
  
  WrapCOASTSS <- function(data) {
    ss <- CalcSumstats(data = data)
    test <- COASTSS(
      anno = ss$sumstats$anno,
      beta = ss$sumstats$beta,
      se = ss$sumstats$se,
      ld = ss$ld,
      maf = ss$sumstats$maf,
      weights = c(1, 2, 3, 4)
    )
    betas <- test@Betas
    return(betas)
  }
  
  withr::local_seed(104) 
  z <- stats::qnorm(0.975)
  
  # Baseline model.
  exp <- c(1, 2, 3, 4)
  data <- DGP(
    beta = exp, 
    method = "none", 
    n = 1e4,
    prop_anno = c(0.25, 0.25, 0.25, 0.25),
    weights = c(1, 1, 1, 1)
  )
  betas <- WrapCOASTSS(data)
  obs <- betas$beta[1:4]
  ses <- betas$se[1:4]
  lower <- obs - z * ses
  upper <- obs + z * ses
  expect_true(all(lower <= exp & exp <= upper))
  
  # Allelic sum model.
  exp <- 1
  data <- DGP(
    beta = exp, 
    method = "sum", 
    n = 1e4,
    prop_anno = c(0.25, 0.25, 0.25, 0.25),
    weights = c(1, 2, 3, 4)
  )
  betas <- WrapCOASTSS(data)
  obs <- betas$beta[betas$test == "sum"]
  ses <- betas$se[betas$test == "sum"]
  lower <- obs - z * ses
  upper <- obs + z * ses
  expect_true(all(lower <= exp & exp <= upper))
  
})

