# Purpose: Data generation for allelic series.
# Updated: 2023-10-04


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

    if (sum(g) == 0) {
      # Randomly add 1 minor allele to avoid MAC = 0.
      draw <- sample.int(n, size = 1)
      g[draw] <- 1
    }

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
#' @param snps Number of SNPs in the gene.
#' @param p_dmv Frequency of deleterious missense variants.
#' @param p_ptv Frequency of protein truncating variants.
#' @return (snps x 1) integer vector.
GenAnno <- function(
  snps,
  p_dmv = 0.33,
  p_ptv = 0.33
) {
  p_bmv <- 1 - sum(c(p_dmv, p_ptv))
  p <- c(p_bmv, p_dmv, p_ptv)
  anno <- stats::rmultinom(snps, size = 1, prob = p)
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
  anno <- GenAnno(snps, p_dmv, p_ptv)
  out <- list(
    anno = anno,
    geno = geno
  )
  return(out)
}


#' Filter Noncausal Variants
#' 
#' Remove a random fraction of variants, which are designated non-causal. 
#' 
#' @param anno (snps x 1) annotation vector.
#' @param geno (n x snps) genotype matrix.
#' @param prop_causal Proportion of variants which are causal.
#' @return List containing the (n x snps) genotype matrix "geno" and the
#' (snps x 1) annotation vector "anno".
FilterGenos <- function(
  anno,
  geno,
  prop_causal = 1.0
) {
  
  n_snp <- length(anno)
  p_noncausal <- 1.0 - prop_causal
  n_noncausal <- round(n_snp * p_noncausal)
  which_drop <- sample(x = n_snp, size = n_noncausal, replace = FALSE)
  
  # Code non-causal variants as "-1". 
  anno[which_drop] <- -1
  geno_filtered <- geno[, anno != -1, drop = FALSE]
  anno_filtered <- anno[anno != -1]
  
  # Output.
  out <- list(
    anno = anno_filtered,
    geno = geno_filtered
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
  N_PC <- 3
  pcs <- lapply(seq_len(N_PC), function(i) {stats::rnorm(n)})
  pcs <- do.call(cbind, pcs)

  # Output.
  out <- cbind(age, pcs)
  out <- scale(out)
  out <- cbind(1, out[, 1], sex, out[, 2:ncol(out)])
  colnames(out) <- c("int", "age", "sex", paste0("pc", seq_len(N_PC)))
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
#' @param prop_causal Proportion of variants which are causal.
#' @param random_signs Randomize signs? FALSE for burden-type genetic
#'   architecture, TRUE for SKAT-type.
#' @param random_var Frailty variance in the case of random signs. Default: 0.
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
  prop_causal = 1.0,
  random_signs = FALSE,
  random_var = 0.0, 
  weights = c(0, 1, 2)
) {

  # Check configuration.
  if ((method != "none") & (length(beta) > 1)) {
    stop("Scalar beta is required for genotype aggregation.")
  }
  if ((method != "none") & random_signs) {
    stop("Random signs are incompatible with aggregated genotypes.")
  }
  
  # Filter genotypes.
  anno_geno <- FilterGenos(
    anno = anno,
    geno = geno,
    prop_causal = prop_causal
  )
  anno <- anno_geno$anno
  geno <- anno_geno$geno

  # Calculate genetic effect.
  # Null phenotype.
  if (all(beta == 0)) {

    eta_g <- 0

  # SKAT-type phenotype.
  } else if(random_signs) {

    n_snps <- length(anno)
    long_beta <- rep(0, n_snps)
    
    # Generate an effect size vector of length == the annotation vector.
    for (i in 0:2) {long_beta[anno == i] <- beta[i + 1]}

    # Sample the sign. 
    signs <- sample(x = c(-1, 1), size = n_snps, replace = TRUE)
    
    # Sample the frailty (positive random effect with expectation 1.0).  
    if (random_var == 0) {
      gamma <- 1
    } else {
      shape <- 1 / random_var
      gamma <- stats::rgamma(n = n_snps, shape = shape, rate = shape)
    }
    
    # Overall genetic effect sizes.
    long_beta <- signs * gamma * long_beta 

    # Convert to indicator genotypes, if required.
    if (indicator) {geno <- 1 * (geno > 0)}
    
    # Overall genetic contribution to the phenotype.
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

    # Overall genetic contribution to the phenotype.
    # If method is not none, beta is a scalar, else a vector.
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

  # Add residual, if required.
  if (include_residual) {
    y <- eta + stats::rnorm(length(eta), sd = reg_param$sd)
  } else {
    y <- eta
  }
  
  # Convert to binary, if required.
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
#' @param anno Annotation vector, if providing genotypes. Should match the
#'   number of columns in geno.
#' @param beta If method = "none", a (3 x 1) coefficient vector for bmvs, dmvs,
#'   and ptvs respectively. If method != "none", a scalar effect size.
#' @param binary Generate binary phenotype? Default: FALSE.
#' @param geno Genotype matrix, if providing genotypes. 
#' @param include_residual Include residual? If FALSE, returns the expected
#'   value. Intended for testing.
#' @param indicator Convert raw counts to indicators? Default: FALSE.
#' @param maf_range Range of minor allele frequencies: c(MIN, MAX).
#' @param method Genotype aggregation method. Default: "none".
#' @param n Sample size.
#' @param p_dmv Frequency of deleterious missense variants. Default of 40% is
#'   based on the frequency of DMVs among rare coding variants in the UK
#'   Biobank.
#' @param p_ptv Frequency of protein truncating variants. Default of 10% is
#'   based on the frequency of PTVs among rare coding variants in the UK
#'   Biobank.
#' @param prop_causal Proportion of variants which are causal. Default: 1.0. 
#' @param random_signs Randomize signs? FALSE for burden-type genetic
#'   architecture, TRUE for SKAT-type.
#' @param random_var Frailty variance in the case of random signs. Default: 0.
#' @param snps Number of SNP in the gene. Default: 100.
#' @param weights Aggregation weights.
#' @return List containing: genotypes, annotations, covariates, phenotypes.
#' @examples 
#' # Generate data.
#' data <- DGP(n = 100)
#' 
#' # View components.
#' table(data$anno)
#' head(data$covar)
#' head(data$geno[, 1:5])
#' hist(data$pheno)
#' @export
DGP <- function(
  anno = NULL,
  beta = c(0, 1, 2),
  binary = FALSE,
  geno = NULL,
  include_residual = TRUE,
  indicator = FALSE,
  maf_range = c(0.005, 0.010),
  method = "none",
  n = 100,
  p_dmv = 0.40,
  p_ptv = 0.10,
  prop_causal = 1.0,
  random_signs = FALSE,
  random_var = 0.0, 
  snps = 100,
  weights = c(1, 2, 3)
) {

  # Generate genotypes and annotations.
  if (is.null(anno) & is.null(geno)) {
    
    # Neither annotations nor genotypes provided.
    anno_geno <- GenGeno(
      n = n,
      snps = snps,
      maf_range = maf_range,
      p_dmv = p_dmv,
      p_ptv = p_ptv
    )
    
  } else if(is.null(anno) & !is.null(geno)) {
    
    # Only genotypes provided.
    geno <- as.matrix(geno)
    n <- nrow(geno)
    snps <- ncol(geno)
    anno <- GenAnno(snps, p_dmv = p_dmv, p_ptv = p_ptv)
    anno_geno <- list(
      anno = anno,
      geno = geno
    )
    
  } else if(!is.null(anno) & is.null(geno)) {
    
    # Only annotations provided.
    snps <- length(anno)
    geno <- GenGenoMat(
      n = n,
      snps = snps,
      maf_range = maf_range
    )
    anno_geno <- list(
      anno = anno,
      geno = geno
    )
    
  } else {
    
    # Annotations and genotypes provided.
    geno <- as.matrix(geno)
    stopifnot(length(anno) == ncol(geno))
    snps <- ncol(geno)
    n <- nrow(geno)
    
    anno_geno <- list(
      anno = anno,
      geno = geno
    )
    
  }

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
    prop_causal = prop_causal,
    random_signs = random_signs,
    random_var = random_var, 
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
