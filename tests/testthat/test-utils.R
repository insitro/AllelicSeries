test_that("Test positive-definite evaluation.", {
  
  # Eigenvalue of zero.
  x <- rbind(c(1, 1), c(1, 1))
  expect_false(isPD(x))  
  
  # Negative eigenvalue (not symmetric).
  x <- rbind(c(1, 1), c(2, 1))
  expect_false(isPD(x))  
  
  # Positive eigenvalues.
  x <- rbind(c(2, 1), c(1, 2))
  expect_true(isPD(x))  
  
})

test_that("Test variant collapsing threshold.", {
  
  # Case 0: Threshold of zero.
  geno <- array(0, dim = c(100, 6))
  anno <- rep(c(1, 2, 3), each = 2)
  colnames(geno) <- paste0("rs", seq_len(6))
  for(i in 1:6) {
    geno[1:(2 * i), i] <- 1
  }
  out <- CollapseGeno(anno = anno, geno = geno, min_mac = 0)
  expect_equal(out$geno, geno)
  
  # Case 1: Threshold of 20.
  # All variants should be collapsed.
  out <- CollapseGeno(anno = anno, geno = geno, min_mac = 20)
  exp <- cbind(
    geno[, 1] + geno[, 2], 
    geno[, 3] + geno[, 4], 
    geno[, 5] + geno[, 6]
  )
  colnames(exp) <- c("agg1", "agg2", "agg3")
  expect_equal(out$geno, exp)
  exp <- c("rs1;rs2", "rs3;rs4", "rs5;rs6")
  expect_equal(out$vars$vars, exp)
  
  # Case 2: Threshold of 7.
  # Variants 1 and 2 are collapsed, variant 3 is "collapsed" with itself.
  # Variants 4, 5, 6 are unchanged.
  out <- CollapseGeno(anno = anno, geno = geno, min_mac = 7)
  exp <- cbind(
    geno[, 1] + geno[, 2], 
    geno[, 4], geno[, 3], geno[, 5], geno[, 6]
  )
  colnames(exp) <- c("agg1", "rs4", "agg2", "rs5", "rs6")
  expect_equal(out$geno, exp)
  exp <- c("rs1;rs2", "rs3", "")
  expect_equal(out$vars$vars, exp)
  
})


test_that("Test variant collapsing with MACs of 0.", {
  
  # Case 0: All zeros provided.
  geno <- array(0, dim = c(100, 3))
  anno <- c(1, 2, 3)
  expect_warning(CollapseGeno(anno = anno, geno = geno))
  
  # Case 1: One category has a MAC of zero.
  geno <- array(0, dim = c(100, 3))
  geno[1, 2] <- geno[1, 3] <- 1
  anno <- c(1, 2, 3)
  out <- suppressWarnings(CollapseGeno(anno = anno, geno = geno, min_mac = 0))
  expect_equal(out$geno, geno[, 2:3])
  
})


test_that("Test variant collapsing with dosage genotypes.", {
  
  # Variants 1 and 2 are collapsed, variant 3 is "collapsed" with itself.
  # Variants 4, 5, 6 are unchanged.
  geno <- array(0, dim = c(100, 6))
  anno <- rep(c(1, 2, 3), each = 2)
  colnames(geno) <- paste0("rs", seq_len(6))
  for(i in 1:6) {
    geno[1:(2 * i), i] <- 0.5
  }
  out <- CollapseGeno(anno = anno, geno = geno, min_mac = 3.5)

  exp <- cbind(
    geno[, 1] + geno[, 2], 
    geno[, 4], geno[, 3], geno[, 5], geno[, 6]
  )
  colnames(exp) <- c("agg1", "rs4", "agg2", "rs5", "rs6")
  expect_equal(out$geno, exp)
  exp <- c("rs1;rs2", "rs3", "")
  expect_equal(out$vars$vars, exp)
  
})

