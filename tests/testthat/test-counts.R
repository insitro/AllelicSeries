test_that("Check variant counter.", {

  # Data.
  anno <- c(1, 2, 3)
  geno <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 2, 0),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 2)
  )
  
  # Counts with minimum MAC = 0.
  counts <- Counts(anno = anno, geno = geno, n_anno = 3, min_mac = 0)
  expect_equal(counts$alleles, c(1, 1 + 2, 1 + 1 + 2))
  expect_equal(counts$variants, c(1, 1, 1))
  expect_equal(counts$carriers, c(1, 2, 3))
  
  # Counts with minimum MAC = 2.
  counts <- Counts(anno = anno, geno = geno, n_anno = 3, min_mac = 2)
  expect_equal(counts$alleles, c(0, 1 + 2, 1 + 1 + 2))
  expect_equal(counts$variants, c(0, 1, 1))
  expect_equal(counts$carriers, c(0, 2, 3))
  
  # Counts with minimum MAC = 100.
  counts <- Counts(anno = anno, geno = geno, n_anno = 3, min_mac = 100)
  expect_equal(counts$alleles, c(0, 0, 0))
  expect_equal(counts$variants, c(0, 0, 0))
  expect_equal(counts$carriers, c(0, 0, 0))
  
  # Check case of multiple variants per annotation category.
  anno <- c(1, 2, 2, 3, 3, 3)
  geno <- diag(6) # Identity matrix.
  counts <- Counts(anno = anno, geno = geno, n_anno = 3, min_mac = 0)
  expect_equal(counts$alleles, c(1, 2, 3))
  expect_equal(counts$variants, c(1, 2, 3))
  expect_equal(counts$carriers, c(1, 2, 3))
  
  # Check case of 4 annotation categories.
  anno <- c(1, 2, 3, 4)
  geno <- rbind(
    c(1, 0, 0, 0),
    c(0, 2, 0, 0),
    c(1, 0, 1, 0),
    c(0, 1, 0, 1)
  )
  counts <- Counts(anno = anno, geno = geno, n_anno = 4)
  expect_equal(counts$alleles, c(2, 3, 1, 1))
  expect_equal(counts$variants, c(1, 1, 1, 1))
  expect_equal(counts$carriers, c(2, 2, 1, 1))
  
  # 4 annotation categories but n_count is set to 3.
  # Counts should be the same as the first 3 values above.
  counts <- Counts(anno = anno, geno = geno, n_anno = 3)
  expect_equal(counts$alleles, c(2, 3, 1))
  expect_equal(counts$variants, c(1, 1, 1))
  expect_equal(counts$carriers, c(2, 2, 1))
  
  # 3 annotation categories but n_count is set to 4.
  anno <- c(1, 2, 3)
  geno <- rbind(
    c(1, 0, 0),
    c(0, 2, 1),
    c(1, 1, 1)
  )
  counts <- Counts(anno = anno, geno = geno, n_anno = 4)
  expect_equal(counts$alleles, c(2, 3, 2, 0))
  expect_equal(counts$variants, c(1, 1, 1, 0))
  expect_equal(counts$carriers, c(2, 2, 2, 0))
})
