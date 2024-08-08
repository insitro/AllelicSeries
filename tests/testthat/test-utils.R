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
