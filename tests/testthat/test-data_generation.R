test_that("Check minor allele count is always positive.", {

  g <- GenGenoMat(n = 100, snps = 1e3)
  mac <- apply(g, 2, sum)
  expect_gte(min(mac), 1)

})


test_that("Check annotation generation.", {

  withr::local_seed(101)
  g <- GenGenoMat(n = 10, snps = 10)
  anno <- GenAnno(ncol(g))
  expect_true(all(anno %in% c(0, 1, 2)))

})


test_that("Check phenotype generation.", {

  anno <- c(0, 1, 2)
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


  # Genotypes for aggregation methods.
  anno <- c(0, 1, 2)
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


test_that("Check data generation with manual phenotypes.", {
  
  # Annotations and genotypes both provided.
  anno <- c(0, 1, 2)
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


test_that("Check ability to vary the proportion of causal variants.", {
  
  anno <- c(0, 1, 2)
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


test_that("Checking ability to generate phenotypes with random genetic effects.", {
  
  anno <- c(0, 1, 2)
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

