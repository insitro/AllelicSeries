test_that("Check genotype generation.", {

  withr::local_seed(101)
  
  # Minor allele count is always positive.
  g <- GenGenoMat(n = 100, snps = 1e3)
  mac <- apply(g, 2, sum)
  expect_gte(min(mac), 1)
  
  # Number of alternate alleles never exceeds 2.
  g <- GenGenoMat(n = 100, snps = 1e3, maf_range = c(0.5, 1.0))
  max_g <- max(g)
  expect_lte(max_g, 2)
  
  # Minor allele frequency is near the expected range.
  # MAF is not expected to fall exactly within the rage due to stochasticity.
  n <- 1e3
  g <- GenGenoMat(n = n, snps = 1e2, maf_range = c(0.001, 0.005))
  maf <- apply(g, 2, mean) / 2
  tol <- 1 / (10 * n) 
  expect_gte(min(maf), 1 / (2 * n) - tol)
  expect_lte(max(maf), 2 * 0.005 + tol)
  
})

# ------------------------------------------------------------------------------

test_that("Check annotation generation.", {

  withr::local_seed(101)
  anno <- GenAnno(1000)
  expect_true(all(anno %in% c(1, 2, 3)))
  
  anno <- GenAnno(1000, prop_anno = c(1, 1, 1, 1))
  expect_true(all(anno %in% c(1, 2, 3, 4)))

})


# ------------------------------------------------------------------------------

test_that("Check phenotype generation.", {
  
  # Test cases without aggregation.
  anno <- c(1, 2, 3)
  geno <- rbind(
    c(0, 0, 0),
    c(2, 0, 0),
    c(0, 2, 0),
    c(0, 0, 2)
  )
  covar <- c(rep(1, 4))
  beta <- c(1, 2, 4)
  reg_param <- list(coef = as.matrix(0))
  weights <- c(1, 1, 1)

  # Counts.
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    indicator = FALSE,
    weights = weights
  )
  exp <- c(0, 2 * 1, 2 * 2, 2 * 4)
  expect_true(all(obs == exp))

  # Indicators.
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    indicator = TRUE,
    weights = weights
  )
  exp <- c(0, 1 * 1, 1 * 2, 1 * 4)
  expect_true(all(obs == exp))

  # Random signs.
  withr::local_seed(101)
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    random_signs = TRUE,
    weights = weights
  )
  exp <- c(0, 2 * -1, 2 * -2, 2 * 4)
  expect_true(all(obs == exp))

  # Binary phenotype.
  withr::local_seed(101)
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    binary = TRUE,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    random_signs = TRUE,
    weights = weights
  )
  exp <- 1 * (c(0, 2 * -1, 2 * -2, 2 * 4) >= 0)
  expect_true(all(obs == exp))
  
  # Test cases with aggregation.
  # Genotypes for aggregation methods.
  anno <- c(1, 2, 3)
  geno <- rbind(
    c(0, 0, 0),
    c(1, 1, 1),
    c(0, 2, 0),
    c(1, 0, 2)
  )
  weights <- c(0, 1, 2)

  # Sum aggregation.
  beta <- 1.0
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    indicator = FALSE,
    method = "sum",
    weights = weights
  )
  exp <- c(0, 1 * 1 + 1 * 2, 2 * 1, 2 * 2)
  expect_equal(obs, exp)

  # Max aggregation.
  beta <- 1.0
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    indicator = FALSE,
    method = "max",
    weights = weights
  )
  exp <- c(0, 2, 2, 4)
  expect_equal(obs, exp)
  
  # Check case of more than 3 annotation categories.
  anno <- c(1, 2, 3, 4)
  geno <- rbind(
    c(0, 0, 0, 0),
    c(2, 0, 0, 0),
    c(0, 2, 0, 0),
    c(0, 0, 2, 0),
    c(0, 0, 0, 2)
  )
  covar <- c(rep(1, 5))
  beta <- c(1, 2, 4, 8)
  reg_param <- list(coef = as.matrix(0))
  weights <- c(1, 1, 1, 1)
  
  # Counts.
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    indicator = FALSE,
    weights = weights
  )
  exp <- c(0, 2 * 1, 2 * 2, 2 * 4, 2 * 8)
  expect_true(all(obs == exp))
  
  # Indicators.
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    indicator = TRUE,
    weights = weights
  )
  exp <- c(0, 1 * 1, 1 * 2, 1 * 4, 1 * 8)
  expect_true(all(obs == exp))
  
  # Genotypes for aggregation.
  geno <- rbind(
    c(0, 0, 0, 0),
    c(1, 0, 0, 0),
    c(0, 1, 0, 0),
    c(1, 0, 1, 0),
    c(0, 1, 0, 1)
  )
  weights <- c(1, 2, 3, 4)
  
  # Sum aggregation.
  beta <- 1.0
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    indicator = FALSE,
    method = "sum",
    weights = weights
  )
  exp <- c(0, 1, 2, 1 + 3, 2 + 4)
  expect_equal(obs, exp)
  
  # Max aggregation.
  obs <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    indicator = FALSE,
    method = "max",
    weights = weights
  )
  exp <- c(0, 1, 2, 3, 4)
  expect_equal(obs, exp)
  
  # Check incompatible settings.
  # Random signs cannot be used with aggregated genotypes.
  expect_error(
    obs <- GenPheno(
      anno = anno,
      beta = beta,
      covar = covar,
      geno = geno,
      reg_param = reg_param,
      random_signs = TRUE,
      method = "max"
    )
  )

})


# ------------------------------------------------------------------------------

test_that("Check data generation with manual phenotypes.", {
  
  # Annotations and genotypes both provided.
  anno <- c(1, 2, 3)
  geno <- rbind(
    c(0, 0, 0),
    c(1, 1, 1),
    c(0, 2, 0),
    c(1, 0, 2)
  )
  
  data <- DGP(anno = anno, geno = geno)
  expect_true(all(data$anno == anno))
  expect_true(all(data$geno == geno))
  
  # Only annotations provided.
  data <- DGP(anno = anno, n = 10, snps = 10)
  expect_true(all(data$anno == anno))
  expect_true(ncol(data$geno) == length(anno))  # snps replaced by length(anno).
  
  # Only genotypes provided.
  data <- DGP(geno = geno, n = 10, snps = 10)
  expect_true(all(data$geno == geno))
  expect_true(length(data$pheno) == nrow(geno))  # n replaced by nrow(geno).
  expect_true(length(data$anno) == ncol(geno))  # snps replaced by ncol(geno).
  
})


# ------------------------------------------------------------------------------

test_that("Check ability to vary the proportion of causal variants.", {
  
  anno <- c(1, 2, 3)
  geno <- rbind(
    c(0, 0, 0),
    c(1, 1, 1),
    c(0, 2, 0),
    c(1, 0, 2)
  )
  covar <- c(rep(1, 4))
  beta <- c(1, 2, 3)
  reg_param <- list(coef = as.matrix(0))
  weights <- c(1, 1, 1)
  
  # Setting the causal proportion to zero results in no genetic component.
  null_pheno <- GenPheno(
    anno = anno,
    beta = beta,
    geno = geno,
    covar = covar,
    reg_param = reg_param,
    include_residual = FALSE,
    prop_causal = 0.0
  )
  expect_true(all(null_pheno == 0))
  
  # Generate phenotypes 2 ways:
  # 1. Filtering genotypes externally.
  # 2. Filtering genotypes internally.
  withr::local_seed(1010)
  filtered_geno <- FilterGenos(anno = anno, geno = geno, prop_causal = 1/3)
  external_pheno <- GenPheno(
    anno = filtered_geno$anno,
    beta = beta,
    geno = filtered_geno$geno,
    covar = covar,
    reg_param = reg_param,
    include_residual = FALSE,
    prop_causal = 1.0
  )
  
  withr::local_seed(1010)
  internal_pheno <- GenPheno(
    anno = anno,
    beta = beta,
    geno = geno,
    covar = covar,
    reg_param = reg_param,
    include_residual = FALSE,
    prop_causal = 1/3
  )
  expect_equal(external_pheno, internal_pheno)
  
})


# ------------------------------------------------------------------------------

test_that(
  "Checking ability to generate phenotypes with random genetic effects.", {
  
  anno <- c(1, 2, 3)
  geno <- rbind(
    c(0, 0, 0),
    c(1, 1, 1),
    c(0, 2, 0),
    c(1, 0, 2)
  )
  covar <- c(rep(1, 4))
  beta <- c(1, 2, 3)
  reg_param <- list(coef = as.matrix(0))
  weights <- c(1, 1, 1)
  
  withr::local_seed(123)
  y_nonrandom <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    random_signs = TRUE,
    random_var = 0.0
  )
  
  withr::local_seed(123)
  y_random <- GenPheno(
    anno = anno,
    beta = beta,
    covar = covar,
    geno = geno,
    reg_param = reg_param,
    include_residual = FALSE,
    random_signs = TRUE,
    random_var = 1.0
  )
  
  expect_gt(stats::var(y_random), stats::var(y_nonrandom))
  
})


# ------------------------------------------------------------------------------

test_that("Check sumstats generation.", {
  
  withr::local_seed(123)
  anno <- c(1, 2, 3)
  n <- 100
  geno <- replicate(3, stats::rbinom(n = n, size = 2, prob = 0.25))
  pheno <- stats::rnorm(n = n)
  covar <- rep(1, n)
  
  # Calculate sumstats manually.
  ld <- cor(geno)
  sumstats <- lapply(seq_len(3), function(i) {
    fit <- stats::lm(pheno ~ geno[, i])
    results <- summary(fit)
    beta <- as.numeric(results$coefficients[, "Estimate"][2])
    se <- as.numeric(results$coefficients[, "Std. Error"][2])
    return(data.frame(beta = beta, se = se))
  })
  sumstats <- do.call(rbind, sumstats)
  
  # Check function outputs.
  data <- list(
    anno = anno, 
    geno = geno, 
    pheno = pheno, 
    covar = covar,
    type = "quantitative"
  )
  obs <- CalcSumstats(data = data)
  
  expect_equal(ld, obs$ld, tolerance = 1e-4)
  expect_equal(sumstats$beta, obs$sumstats$beta, tolerance = 1e-4)
  expect_equal(sumstats$se, obs$sumstats$se, tolerance = 1e-4)
  
})


test_that("Test addition of intercept during sumstat calculation.", {
  
  withr::local_seed(123)
  data <- DGP()
  
  # Results when an intercept is included manually.
  with_int <- CalcSumstats(data = data)
  
  # Results without intercept.
  covar <- data$covar
  covar <- covar[, 2:6]
  expect_warning({
    without_int <- CalcSumstats(
      anno = data$anno,
      geno = data$geno,
      covar = covar,
      pheno = data$pheno
    )
  })
  
  expect_equal(
    with_int$sumstats, 
    without_int$sumstats
  )
  
})


# ------------------------------------------------------------------------------

test_that("Test sumstat generation when covariates are omitted.", {
  
  withr::local_seed(123)
  data <- DGP()
  
  # Method 1: No covariates provided.
  sumstats1 <- CalcSumstats(
    anno = data$anno, 
    geno = data$geno, 
    pheno = data$pheno
  )
  
  # Method 2: Provide intercept for covariates.
  sumstats2 <- CalcSumstats(
    anno = data$anno, 
    covar = rep(1, length(data$pheno)),
    geno = data$geno, 
    pheno = data$pheno
  )
  
  expect_equal(sumstats1$sumstats, sumstats2$sumstats)
  
})


# ------------------------------------------------------------------------------

test_that(
  "Expect error when the number of annotation categories is inconsistent.", {
    
    # Fewer weights than annotation categories.
    expect_error(
      data <- DGP(
        prop_anno = c(1, 1, 1, 1),
        weights = c(1, 2, 3)
      )
    )
    
    # Fewer betas than annotation categories.
    expect_error(
      data <- DGP(
        beta = c(1, 2, 3),
        prop_anno = c(1, 1, 1, 1),
        weights = c(1, 2, 3, 4)
      )
    )
    
    # Expect no error because parameters are consistent.
    expect_error(
      data <- DGP(
        beta = c(1, 2, 3, 4),
        prop_anno = c(1, 1, 1, 1),
        weights = c(1, 2, 3, 4)
      ), NA)
    
    # Expect no error because aggregation is applied.
    expect_error(
      data <- DGP(
        beta = 1,
        prop_anno = c(1, 1, 1, 1),
        method = "sum",
        weights = c(1, 2, 3, 4)
      ), NA)
    
  })
