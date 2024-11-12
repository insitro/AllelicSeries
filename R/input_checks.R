# Purpose: Input checks.
# Updated: 2024-11-11

#' Relevel Annotations
#' 
#' @param anno (snps x 1) annotation vector.
#' @return (snps x 1) annotation vector.
#' @noRd
RelevelAnno <- function(anno) {
  int_anno <- as.integer(anno)
  if (any(is.na(int_anno)) | !all(anno == int_anno)) {
    stop("The annotation vector should use integer labels of 1:length(weights).")
  }
  if (0 %in% anno) {
    warning(
      "Zero is present in the annotation vector. Annotations will be reindexed
      to start at one."
    )
    anno <- anno + 1
  }
  return(anno)
}


#' Check Inputs
#'
#' @param anno (snps x 1) annotation vector.
#' @param covar (n x p) covariate matrix.
#' @param geno (n x snps) genotype matrix.
#' @param is_pheno_binary Is the phenotype binary?
#' @param pheno (n x 1) phenotype vector.
#' @param weights (L x 1) annotation category weights.
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
  
  # Check annotation values.
  n_anno <- length(weights)
  if (!all(anno %in% seq_len(n_anno))) {
    stop("The annotations should take values in 1:length(weights).")
  }
  if (!all(seq_len(n_anno) %in% anno)) {
    warning("Annotation categories are present to which no variants are assigned.")
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
#' @param lambda Genomic inflation factor.
#' @param ld (snps x snps) matrix of correlations among the genetic variants.
#'   Although ideally provided, an identity matrix is assumed if not.
#' @param weights (L x 1) annotation category weights.
#' @param is_skat Logical, is the check for the SKAT test?
#' @param maf (snps x 1) vector of minor allele frequencies. Although ideally
#'   provided, defaults to the zero vector.
#' @return Logical indicating whether the matrix was positive definite.
CheckInputsSS <- function(
    anno,
    beta,
    se,
    lambda,
    ld,
    weights,
    is_skat = FALSE,
    maf = NULL
) {
  
  # Check annotation values.
  n_anno <- length(weights)
  if (!all(anno %in% seq_len(n_anno))) {
    stop("The annotations should take values in 1:length(weights).")
  }
  if (!all(seq_len(n_anno) %in% anno)) {
    warning("Annotation categories are present to which no variants are assigned.")
  }
  
  # Check inflation factor is >= 1.
  if (any(lambda < 1)) {
    msg <- paste0(
      "The inflation factor labmda should be >= 1. ",
      "Values < 1 will be reset to 1."
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
  out <- TRUE
  if (!is.null(ld)) {
    is_pd <- isPD(ld)
    if (!is_pd) {
      msg <- paste0(
        "LD has very small or negative eigenvalues. ",
        "Epsilon will be added to the diagonal."
      )
      warning(msg)
      out <- FALSE
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
  if (is_skat & !is.null(maf)) {any_na <- any_na | any(is.na(maf))}
  
  if (any_na) {
    stop("Input data should not contain missing values.")
  }
  
  # Output indicates whether the LD matrix is PSD.
  return(out)
}
