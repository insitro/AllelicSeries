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
  
  # Counts.
  counts <- Counts(anno = anno, geno = geno, min_mac = 0)
  expect_equal(counts$alleles, c(1, 1 + 2, 1 + 1 + 2))
  expect_equal(counts$variants, c(1, 1, 1))
  expect_equal(counts$carriers, c(1, 2, 3))
  
  # Counts with minimum MAC.
  counts <- Counts(anno = anno, geno = geno, min_mac = 2)
  expect_equal(counts$alleles, c(0, 1 + 2, 1 + 1 + 2))
  expect_equal(counts$variants, c(0, 1, 1))
  expect_equal(counts$carriers, c(0, 2, 3))
  
  # Check case of multiple varaints per annotation category.
  anno <- c(0, 1, 1, 2, 2, 2)
  geno <- diag(6) # Identity matrix.
  counts <- Counts(anno = anno, geno = geno, min_mac = 0)
  expect_equal(counts$alleles, c(1, 2, 3))
  expect_equal(counts$variants, c(1, 2, 3))
  expect_equal(counts$carriers, c(1, 2, 3))

})
