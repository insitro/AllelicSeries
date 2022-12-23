# Purpose: Input checks.
# Updated: 2022-09-22

#' Check Inputs
#'
#' @param anno (snps x 1) annotation vector.
#' @param covar (n x p) covariate matrix.
#' @param geno (n x snps) genotype matrix.
#' @param is_pheno_binary Is the phenotype binary?
#' @param pheno (n x 1) phenotype vector.
#' @param weights (3 x 1) annotation category weights.
#' @return None.
CheckInputs <- function(
    anno,
    covar,
    geno,
    is_pheno_binary,
    pheno,
    weights
) {

  # Check annotation length.
  if (length(anno) != ncol(geno)) {
    stop("Length of the annotation vector should match number of columns in
         the genotype matrix.")
  }

  # Check for missing values.
  any_na <- any(is.na(anno)) |
    any(is.na(covar)) |
    any(is.na(geno)) |
    any(is.na(pheno))

  if (any_na) {
    stop("Input data should not contain missing values.")
  }

  # Check dimensional consistency.
  equal_dims <- all.equal(nrow(covar), nrow(geno), length(pheno))
  if (!equal_dims) {
    stop("Input dimensions are inconsistent.")
  }

  # Check for non-zero genotypes and phenotypes.
  if (all(geno == 0)) {
    stop("No non-zero genotypes are present.")
  }
  if (all(pheno == 0)){
    stop("No non-zero phenotypes are present.")
  }

    # Check weights.
    if (all(weights == 0)) {
      stop("At least 1 annotation category requires non-zero weights.")
    }
  if (any(weights < 0)) {
    stop("Annotation category weights should be non-negative.")
  }

  # Check binary phenotype.
  unique_pheno_values <- sort(unique(pheno))
  if (is_pheno_binary) {
    improper <- any(
      length(unique_pheno_values) != 2,
      !all(unique_pheno_values %in% c(0, 1))
    )
    if (improper) {
      stop("A binary phenotype must have exactly 2 values: {0, 1}.")
    }
  }
}
