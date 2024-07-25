test_that("Check indicator aggregation.", {

  # Method = "none".
  weights <- c(1, 1, 1)
  anno <- c(0, 1, 2)
  geno <- rbind(c(0, 0, 0), c(3, 0, 0), c(0, 2, 0), c(0, 1, 1))
  obs <- Aggregator(
    anno, geno, indicator = TRUE, method = "none", weights = weights)
  exp <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 1, 1))
  expect_true(all(obs == exp))

  # Method = "sum".
  weights <- c(0, 1, 2)
  obs <- Aggregator(
    anno, geno, indicator = TRUE, method = "sum", weights = weights)
  exp <- c(0, 0, 1, 1 + 2)
  expect_true(all(obs == exp))

  # Method = "max".
  weights <- c(0, 1, 2)
  obs <- Aggregator(
    anno, geno, indicator = TRUE, method = "max", weights = weights)
  exp <- c(0, 0, 1, 2)
  expect_true(all(obs == exp))

})


test_that("Check count aggregation.", {

  # Method = "none".
  weights <- c(1, 1, 1)
  anno <- rep(c(0, 1, 2), each = 2)
  geno <- rbind(
    c(0, 0, 0, 0, 0, 0),
    c(0, 1, 0, 2, 0, 1),
    c(0, 0, 1, 1, 2, 1)
  )
  obs <- Aggregator(
    anno, geno, method = "none", weights = weights)
  exp <- rbind(c(0, 0, 0), c(1, 2, 1), c(0, 2, 3))
  expect_true(all(obs == exp))

  # Method = "sum".
  weights <- c(0, 1, 2)
  obs <- Aggregator(anno, geno, method = "sum", weights = weights)
  exp <- c(0, 2 + 2, 2 + 3 * 2)
  expect_true(all(obs == exp))

  # Method = "max".
  weights <- c(0, 1, 2)
  obs <- Aggregator(anno, geno, method = "max", weights = weights)
  exp <- c(0, 2, 3 * 2)
  expect_true(all(obs == exp))

})


test_that("Check that omnibus test runs.", {

  anno <- c(0, 1, 2)
  geno <- rbind(
    c(0, 0, 0),
    c(2, 0, 0),
    c(0, 2, 0),
    c(0, 0, 2),
    c(1, 1, 1)
  )
  pheno <- c(0, 1, 2, 3, 4)

  # No error expected.
  expect_error(COAST(anno, geno, pheno), NA)

})


test_that("Case of all zero weights.", {

  anno <- c(0, 1, 2)
  geno <- rbind(
    c(0, 0, 0),
    c(0, 0, 1),
    c(1, 0, 0)
  )
  pheno <- c(0, 1, 2)
  weights <- c(0, 0, 0)

  expect_warning(Aggregator(anno, geno, method = "none", weights = weights))
  expect_error(ASBT(anno, geno, pheno, weights = weights))
  expect_error(ASKAT(anno, geno, pheno, weights = weights))
  expect_error(COAST(anno, geno, pheno, weights = weights))

})


test_that("Case of all zero genotypes", {

  anno <- c(0, 1, 2)
  geno <- rbind(
    c(0, 0, 0),
    c(0, 0, 0),
    c(0, 0, 0)
  )
  pheno <- c(0, 1, 2)

  expect_warning(Aggregator(anno, geno, method = "none"))
  expect_error(ASBT(anno, geno, pheno))
  expect_error(ASKAT(anno, geno, pheno))
  expect_error(COAST(anno, geno, pheno))

})


test_that("Validate input checks.", {

  anno <- c(0, 1, 2)
  covar <- c(1, 1, 1)
  geno <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  pheno <- c(0, 1, 2)
  weights <- c(0, 1, 2)

  # Error expected: all genotypes zero.
  expect_error(
    CheckInputs(
      anno = anno,
      covar = covar,
      geno = 0 * geno,
      pheno = pheno,
      weights = weights
    )
  )

  # Error expected: all phenotypes zero.
  expect_error(
    CheckInputs(
      anno = anno,
      covar = covar,
      geno = geno,
      pheno = 0 * pheno,
      weights = weights
    )
  )

  # Error expected: phenotype is not binary.
  expect_error(
    CheckInputs(
      anno = anno,
      covar = covar,
      geno = geno,
      is_pheno_binary = TRUE,
      pheno = pheno,
      weights = weights
    )
  )

  # Error expected: negative weights.
  expect_error(
    CheckInputs(
      anno = anno,
      covar = covar,
      geno = geno,
      is_pheno_binary = TRUE,
      pheno = 1 * (pheno > 0),
      weights = c(-1, 0, 1)
    )
  )

})


test_that("Overall check of omnibus test.", {

  anno <- c(0, 1, 2)
  geno <- rbind(
    c(0, 0, 0),
    c(0, 0, 1),
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  n <- nrow(geno)
  pheno <- c(-1, 1, 0, 0, 2)

  # Without standard SKAT-O all or PTV.
  p_omni <- COAST(
    anno = anno,
    geno = geno,
    pheno = pheno,
    include_orig_skato_all = FALSE,
    include_orig_skato_ptv = FALSE
  )
  pvals <- p_omni@Pvals
  p_omni <- as.numeric(pvals$pval[pvals$test == "omni"])
  expect_equal(length(pvals$pval), 8)
  expect_equal(p_omni, 0.0, tolerance = 0.005)

  # With SKAT-O all.
  p_omni <- COAST(
    anno = anno,
    geno = geno,
    pheno = pheno,
    include_orig_skato_all = TRUE,
    include_orig_skato_ptv = FALSE
  )
  pvals <- p_omni@Pvals
  p_omni <- as.numeric(pvals$pval[pvals$test == "omni"])
  expect_equal(length(pvals$pval), 8 + 1)
  expect_equal(p_omni, 0.0, tolerance = 0.005)

  # With SKAT-O PTV.
  p_omni <- COAST(
    anno = anno,
    geno = geno,
    pheno = pheno,
    include_orig_skato_all = FALSE,
    include_orig_skato_ptv = TRUE
  )
  pvals <- p_omni@Pvals
  p_omni <- as.numeric(pvals$pval[pvals$test == "omni"])
  expect_equal(length(pvals$pval), 8 + 1)
  expect_equal(p_omni, 0.0, tolerance = 0.005)

  # With both SKAT-O all and SKAT-O PTV.
  p_omni <- COAST(
    anno = anno,
    geno = geno,
    pheno = pheno,
    include_orig_skato_all = TRUE,
    include_orig_skato_ptv = TRUE
  )
  pvals <- p_omni@Pvals
  p_omni <- as.numeric(pvals$pval[pvals$test == "omni"])
  expect_equal(length(pvals$pval), 8 + 2)
  expect_equal(p_omni, 0.0, tolerance = 0.005)

})


test_that("Check application to common variants.", {
  
  anno <- c(0, 1, 2)
  geno <- rbind(
    c(1, 2, 1),
    c(1, 1, 2),
    c(1, 2, 2),
    c(0, 1, 1),
    c(0, 2, 2)
  )
  n <- nrow(geno)
  pheno <- c(-1, 1, 0, 3, 2)
  
  expect_warning({
    COAST(
      anno = anno,
      geno = geno,
      pheno = pheno
    )
  })
  
})
