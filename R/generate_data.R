# Purpose: Data generation for allelic series.
# Updated: 2022-08-02


#' Generate Genotype Matrix
#'
#' @param n Sample size.
#' @param snps Number of SNP in the gene.
#' @param maf_range Range of minor allele frequencies: c(MIN, MAX).
#' @return (n x snps) numeric matrix.
GenGenoMat <- function(
  n,
  snps,
  maf_range = c(0.005, 0.010)
) {

  geno <- lapply(seq_len(snps), function(i) {

    maf <- stats::runif(1, maf_range[1], maf_range[2])
    g <- stats::rbinom(n, 2, maf)

    # Randomly add 1 minor allele to avoid MAC = 0.
    draw <- sample.int(n, size = 1)
    g[draw] <- g[draw] + 1

    return(g)
  })
  geno <- do.call(cbind, geno)

  return(geno)
}


#' Generate Genotype Annotations
#'
#' Returns a vector of length = the number of columns (SNPs) in the genotype
#' matrix. Each SNP is classified as a benign missense variant (0), a
#' deleterious missense variant (1), or a protein truncating variant (2).
#'
#' @param mat Genotype matrix.
#' @param p_dmv Frequency of deleterious missense variants.
#' @param p_ptv Frequency of protein truncating variants.
#' @return (snps x 1) integer vector.
GenAnno <- function(
  mat,
  p_dmv = 0.33,
  p_ptv = 0.33
) {
  p_bmv <- 1 - sum(c(p_dmv, p_ptv))
  p <- c(p_bmv, p_dmv, p_ptv)
  n_snps <- ncol(mat)

  anno <- stats::rmultinom(n_snps, size = 1, prob = p)
  anno <- apply(anno, 2, which.max) - 1
  return(anno)
}


#' Generate Genotypes
#'
#' @param n Sample size.
#' @param snps Number of SNP in the gene.
#' @param maf_range Range of minor allele frequencies: c(MIN, MAX).
#' @param p_dmv Frequency of deleterious missense variants.
#' @param p_ptv Frequency of protein truncating variants.
#' @return List containing the (n x snps) genotype matrix "geno" and the
#' (snps x 1) annotation vector "anno".
GenGeno <- function(
  n,
  snps,
  maf_range = c(0.005, 0.010),
  p_dmv = 0.33,
  p_ptv = 0.33
) {
  geno <- GenGenoMat(n = n, snps = snps, maf_range = maf_range)
  anno <- GenAnno(geno, p_dmv, p_ptv)
  out <- list(
    geno = geno,
    anno = anno
  )
  return(out)
}


#' Generate Covariates
#'
#' Generate an (n x 6) covariate matrix with columns representing an intercept,
#' age, sex, and 3 genetic PCs. Because these simulations address rare variant
#' analysis, correlation between genotypes and the genetic PCs (based on
#' common variants) is unnecessary.
#'
#' @param n Sample size.
#' @return (n x 6) numeric matrix.
GenCovar <- function(n) {

  # Covariate representing (standardized) age.
  age <- stats::rnorm(n)

  # Covariate representing sex.
  sex <- stats::rbinom(n, size = 1, prob = 0.5)

  # 3 covariates representing PCs.
  pcs <- lapply(seq_len(3), function(i) {stats::rnorm(n)})
  pcs <- do.call(cbind, pcs)

  # Output.
  out <- cbind(age, sex, pcs)
  out <- scale(out)
  out <- cbind(1, out)
  colnames(out) <- c("int", "age", "sex", "pc1", "pc2", "pc3")
  return(out)
}


#' Calculate Regression Parameters
#'
#' Calculate phenotypic regression coefficients and the residual variation
#' based on proportion of variation explained (PVE) by each factor. Note that
#' the proportion of variation explained by genotype is required, but genetic
#' effects are not generated here.
#'
#' @param pve_age PVE by age.
#' @param pve_pcs PVE by PCs (collectively).
#' @param pve_sex PVE by sex.
#' @return List containing the (5 x 1) regression coefficient vector "coef" and
#' the residual standard deviation "sd".
CalcRegParam <- function(
  pve_age = 0.10,
  pve_pcs = 0.20,
  pve_sex = 0.10
) {

  # Residual PVE.
  pve_res <- 1 - sum(c(pve_age, pve_pcs, pve_sex))

  # Draw regression coefficients.
  b_int <- 0
  b_age <- stats::rnorm(1, sd = sqrt(pve_age))
  b_sex <- stats::rnorm(1, sd = sqrt(pve_sex))
  b_pcs <- stats::rnorm(3, sd = sqrt(pve_pcs / 3))
  beta <- c(b_int, b_age, b_sex, b_pcs)
  names(beta) <- c("int", "age", "sex", "pc1", "pc2", "pc3")

  # Output.
  out <- list(
    coef = beta,
    sd = sqrt(pve_res)
  )
  return(out)
}


#' Generate Phenotypes
#'
#' @param anno (snps x 1) annotation vector.
#' @param beta (3 x 1) coefficient vector for bmvs, dmvs, and ptvs respectively.
#' @param covar Covariate matrix.
#' @param geno (n x snps) genotype matrix.
#' @param reg_param Regression parameters.
#' @param binary Generate binary phenotype? Default: FALSE.
#' @param include_residual Include residual? If FALSE, returns the expected
#'   value. Intended for testing.
#' @param indicator Convert raw counts to indicators? Default: FALSE.
#' @param method Genotype aggregation method. Default: "none".
#' @param random_signs Randomize signs? FALSE for burden-type genetic
#'   architecture, TRUE for SKAT-type.
#' @param weights Aggregation weights.
#' @return (n x 1) numeric vector.
GenPheno <- function(
  anno,
  beta,
  covar,
  geno,
  reg_param,
  binary = FALSE,
  include_residual = TRUE,
  indicator = FALSE,
  method = "none",
  random_signs = FALSE,
  weights = c(0, 1, 2)
) {

  # Check configuration.
  if ((method != "none") & (length(beta) > 1)) {
    stop("Scalar beta is required for genotype aggregation.")
  }
  if ((method != "none") & random_signs) {
    stop("Random signs are incompatible with aggregated genotypes.")
  }

  # Calculate genetic effect.
  # Null phenotype.
  if (all(beta == 0)) {

    eta_g <- 0

  # SKAT-type phenotype.
  } else if(random_signs) {

    long_beta <- rep(0, length(anno))
    for (i in 0:2) {long_beta[anno == i] <- beta[i + 1]}

    signs <- sample(c(-1, 1), length(anno), TRUE)
    long_beta <- long_beta * signs

    if (indicator) {geno <- 1 * (geno > 0)}

    eta_g <- as.numeric(geno %*% long_beta)

  # Burden-type phenotypes.
  } else {

    agg_geno <- Aggregator(
      anno = anno,
      geno = geno,
      drop_empty = FALSE,
      indicator = indicator,
      method = method,
      weights = weights
    )

    if (method != "none") {
      eta_g <- as.numeric(agg_geno * beta)
    } else {
      eta_g <- as.numeric(agg_geno %*% beta)
    }

  }

  # Covariate effect.
  eta_x <- as.numeric(covar %*% reg_param$coef)

  # Linear predictor.
  eta <- eta_g + eta_x

  # Final phenotype.
  if (include_residual) {
    y <- eta + stats::rnorm(length(eta), sd = reg_param$sd)
  } else {
    y <- eta
  }
  if (binary) {
    y <- 1 * (y >= 0)
  }

  return(y)
}


# -----------------------------------------------------------------------------


#' Data Generating Process
#' 
#' Generate a data set consisting of: \itemize{
#' \item{"anno"}{A SNP-length annotation vector.}
#' \item{"covar"}{A subject by 6 covariate matrix.}
#' \item{"geno"}{A subject by SNP genotype matrix.}
#' \item{"pheno"}{A subject-length phenotype vector.}
#' }
#'
#' @param n Sample size.
#' @param snps Number of SNP in the gene.
#' @param beta If method = "none", a (3 x 1) coefficient vector for bmvs, dmvs,
#'   and ptvs respectively. If method != "none", a scalar effect size.
#' @param binary Generate binary phenotype? Default: FALSE.
#' @param include_residual Include residual? If FALSE, returns the expected
#'   value. Intended for testing.
#' @param indicator Convert raw counts to indicators? Default: FALSE.
#' @param maf_range Range of minor allele frequencies: c(MIN, MAX).
#' @param method Genotype aggregation method. Default: "none".
#' @param random_signs Randomize signs? FALSE for burden-type genetic
#'   architecture, TRUE for SKAT-type.
#' @param weights Aggregation weights.
#' @return List containing: genotypes, annotations, covariates, phenotypes.
#' @examples 
#' # Generate data.
#' data <- DGP(n = 100, snps = 20)
#' 
#' # View components.
#' table(data$anno)
#' head(data$covar)
#' head(data$geno[, 1:5])
#' hist(data$pheno)
#' @export
DGP <- function(
  n,
  snps,
  beta = c(0, 1, 2),
  binary = FALSE,
  include_residual = TRUE,
  indicator = FALSE,
  maf_range = c(0.005, 0.010),
  method = "none",
  random_signs = FALSE,
  weights = c(0, 1, 2)
) {

  # Generate genotypes.
  anno_geno <- GenGeno(n = n, snps = snps, maf_range = maf_range)

  # Generate covariates.
  covar <- GenCovar(n)

  # Regression parameters.
  reg_param <- CalcRegParam()

  # Generate phenotypes.
  pheno <- GenPheno(
    anno = anno_geno$anno,
    beta = beta,
    covar = covar,
    geno = anno_geno$geno,
    reg_param = reg_param,
    binary = binary,
    include_residual = include_residual,
    indicator = indicator,
    method = method,
    random_signs = random_signs,
    weights = weights
  )

  # Output.
  out <- list(
    anno = anno_geno$anno,
    covar = covar,
    geno = anno_geno$geno,
    pheno = pheno
  )
  return(out)
}
