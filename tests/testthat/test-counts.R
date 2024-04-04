test_that("Check variant counter.", {

  # Data.
  anno <- c(0, 1, 2)
  geno <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 2, 0),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 2)
  )
  
  # No MAC filter: count alleles.
  obs <- CountAlleles(anno, geno, count_carriers = FALSE, min_mac = 0)
  exp <- c(n_bmv = 1, n_dmv = 3, n_ptv = 4)
  expect_equal(obs, exp)
  
  # No MAC filter: count carriers.
  obs <- CountAlleles(anno, geno, count_carriers = TRUE, min_mac = 0)
  exp <- c(n_bmv = 1, n_dmv = 2, n_ptv = 3)
  expect_equal(obs, exp)
  
  # MAC filter of 2: count alleles.
  obs <- CountAlleles(anno, geno, count_carriers = FALSE, min_mac = 2)
  exp <- c(n_bmv = 0, n_dmv = 3, n_ptv = 4)
  expect_equal(obs, exp)
  
  # MAC filter of 2: count carriers.
  obs <- CountAlleles(anno, geno, count_carriers = TRUE, min_mac = 2)
  exp <- c(n_bmv = 0, n_dmv = 2, n_ptv = 3)
  expect_equal(obs, exp)
  
  # MAC filter of 100: count alleles.
  obs <- CountAlleles(anno, geno, count_carriers = FALSE, min_mac = 100)
  exp <- c(n_bmv = 0, n_dmv = 0, n_ptv = 0)
  expect_equal(obs, exp)
  
  # MAC filter of 100: count carriers.
  obs <- CountAlleles(anno, geno, count_carriers = TRUE, min_mac = 100)
  exp <- c(n_bmv = 0, n_dmv = 0, n_ptv = 0)
  expect_equal(obs, exp)

})
