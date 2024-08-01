# Purpose: Input checks.
# Updated: 2024-07-31

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


#' Input Checks for Summary Statistics
#'
#' @param anno (snps x 1) annotation vector with values in c(0, 1, 2).
#' @param beta (snps x 1) vector of effect sizes for the coding genetic variants
#'   within a gene.
#' @param se (snps x 1) vector of standard errors for the effect sizes.
#' @param maf (snps x 1) vector of minor allele frequencies. Although ideally
#'   provided, defaults to the zero vector.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Although ideally provided, an identity matrix is assumed if not.
#' @return None
CheckInputsSS <- function(
    anno,
    beta,
    se,
    maf,
    ld
) {
  
  # Raise warnings if MAF is omitted.
  if (is.null(maf)) {
    msg <- paste0(
      "If MAF is not provided, a zero vector is assumed. ",
      "This may not be accurate in cases where the MAF is appreciably > 0."
    )
    warning(msg)
  }
  
  # Raise warnings if LD is omitted.
  if (is.null(ld)) {
    msg <- paste0(
      "If LD is not provided, an identity matrix is assumed. ",
      "This may not be accurate in cases where the LD is appreciable."
    )
    warning(msg)
  }
  
  # Check that the LD matrix is SPD.
  if (!is.null(ld)) {
    lambda <- base::eigen(x = ld, symmetric = TRUE, only.values = TRUE)
    if (min(lambda$values) <= 1e-8) {
      msg <- paste0(
        "LD has very small or negative eigenvalues. ",
        "Consider adding a small positive constant to the diagonal."
      )
      warning(msg)
    }
  }
  
  # Check dimensions.
  n_anno <- length(anno)
  n_beta <- length(beta)
  n_se <- length(se)
  if (is.null(ld)) {n_ld <- n_anno} else {n_ld <- nrow(ld)}
  if (is.null(maf)) {n_maf <- n_anno} else {n_maf <- length(maf)}
  
  if (!all.equal(n_anno, n_beta, n_se, n_maf)) {
    msg <- paste0(
      "anno, beta, maf, and se should all have the same length. ",
      "ld should have dimensions of length(beta) x length(beta)."
    )
    stop(msg)
  }
  
  # Check for missing values.
  any_na <- any(is.na(anno)) |
    any(is.na(beta)) |
    any(is.na(se)) 
  if (!is.null(ld)) {any_na <- any_na | any(is.na(ld))}
  if (!is.null(maf)) {any_na <- any_na | any(is.na(maf))}
  
  if (any_na) {
    stop("Input data should not contain missing values.")
  }
  
}
