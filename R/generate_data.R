# Purpose: Data generation for allelic series.
# Updated: 2024-11-11


#' Generate Genotype Matrix
#' 
#' Generate genotypes for `n` subject at `snps` variants in linkage 
#' equilibrium. Genotypes are generated such that the MAC is always >= 1.
#'
#' @param n Sample size.
#' @param snps Number of SNP in the gene.
#' @param maf_range Range of minor allele frequencies: c(MIN, MAX).
#' @return (n x snps) numeric matrix.
GenGenoMat <- function(
  n,
  snps,
  maf_range = c(0.001, 0.005)
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
#' matrix. Each SNP is categorized into one of L categories, where L is 
#' determined by the length of `prop_anno`. 
#'
#' @param snps Number of SNPs in the gene.
#' @param prop_anno Proportions of annotations in each category. Length should
#'   equal the number of annotation categories. Default of c(0.5, 0.4, 0.1) is
#'   based on the approximate empirical frequencies of BMVs, DMVs, and PTVs. 
#' @return (snps x 1) integer vector.
GenAnno <- function(
  snps,
  prop_anno = c(0.5, 0.4, 0.1)
) {
  
  # Ensure proportions sum to 1.
  prop_anno <- prop_anno / sum(prop_anno)
  
  # Generate multinomial annotations.
  anno <- stats::rmultinom(snps, size = 1, prob = prop_anno)
  anno <- apply(anno, 2, which.max)
  return(anno)
}


#' Generate Genotypes
#' 
#' Generates genotypes in linkage equilibrium with accompanying annotations.
#'
#' @param n Sample size.
#' @param snps Number of SNP in the gene.
#' @param maf_range Range of minor allele frequencies: c(MIN, MAX).
#' @param prop_anno Proportions of annotations in each category. Length should
#'   equal the number of annotation categories. Default of c(0.5, 0.4, 0.1) is
#'   based on the approximate empirical frequencies of BMVs, DMVs, and PTVs.
#' @return List containing the (n x snps) genotype matrix "geno" and the
#' (snps x 1) annotation vector "anno".
GenGeno <- function(
  n,
  snps,
  maf_range = c(0.001, 0.005),
  prop_anno = c(0.5, 0.4, 0.1)
) {
  geno <- GenGenoMat(n = n, snps = snps, maf_range = maf_range)
  anno <- GenAnno(snps = snps, prop_anno = prop_anno)
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
#' based on proportion of variation explained (PVE) by each factor.
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
#' Simulate a phenotype based on annotations, covariates, and genotypes.
#' 
#' @section Phenotype generation:
#' * To generate phenotypes from the baseline model, set `method` to "none" and
#'   provide a vector `beta` of length equal to the number of annotation 
#'   categories specifying the effect sizes of each.
#' * To generate phenotypes from the allelic series burden models, set `method` 
#'   to "max" or "sum" and provide a scalar `beta.`
#' * To generate phenotypes from the allelic series SKAT model, set `method` to
#'   "none", set `random_signs` to true, and provide a vector `beta` of length
#'   equal to the number of annotation categories.
#'
#' @param anno (snps x 1) annotation vector.
#' @param beta If method = "none", a (L x 1) coefficient with effect sizes for
#'   each annotation category. By default, there are L = 3 annotation categories
#'   corresponding to BMVs, DMVs, and PTVs. If method != "none", a scalar
#'   effect size for the allelic series burden score.
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
#' @param weights Annotation category weights used for aggregation if 
#'   method != "none". 
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
  weights = c(1, 1, 1)
) {

  # Check configuration.
  if ((method != "none") & (length(beta) > 1)) {
    stop("Scalar beta is required for genotype aggregation.")
  }
  if ((method != "none") & random_signs) {
    stop("Random signs are incompatible with aggregated genotypes.")
  }
  if ((method == "none") & (length(beta) != length(weights))) {
    stop("When no aggregation is applied, length(beta) must equal length(weights).")
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
    n_anno <- length(weights)
    long_beta <- rep(0, n_snps)
    
    # Generate an effect size vector of length == the annotation vector.
    for (i in seq_len(n_anno)) {long_beta[anno == i] <- beta[i]}

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
#' Generate a data set consisting of: 
#' * `anno`: (snps x 1) annotation vector.
#' * `covar`: (subjects x 6) covariate matrix.
#' * `geno`: (subjects x snps) genotype matrix.
#' * `pheno`: (subjects x 1) phenotype vector.
#' * `type`: Either "binary" or "quantitative".
#'
#' @param anno Annotation vector, if providing genotypes. Should match the
#'   number of columns in geno.
#' @param beta If method = "none", a (L x 1) coefficient with effect sizes for
#'   each annotation category. By default, there are L = 3 annotation categories
#'   corresponding to BMVs, DMVs, and PTVs. If method != "none", a scalar
#'   effect size for the allelic series burden score.
#' @param binary Generate binary phenotype? Default: FALSE.
#' @param geno Genotype matrix, if providing genotypes. 
#' @param include_residual Include residual? If FALSE, returns the expected
#'   value. Intended for testing.
#' @param indicator Convert raw counts to indicators? Default: FALSE.
#' @param maf_range Range of minor allele frequencies: c(MIN, MAX).
#' @param method Genotype aggregation method. Default: "none".
#' @param n Sample size.
#' @param prop_anno Proportions of annotations in each category. Length should
#'   equal the number of annotation categories. Default of c(0.5, 0.4, 0.1) is
#'   based on the approximate empirical frequencies of BMVs, DMVs, and PTVs. 
#' @param prop_causal Proportion of variants which are causal. Default: 1.0. 
#' @param random_signs Randomize signs? FALSE for burden-type genetic
#'   architecture, TRUE for SKAT-type.
#' @param random_var Frailty variance in the case of random signs. Default: 0.
#' @param snps Number of SNP in the gene. Default: 100.
#' @param weights Annotation category weights. Length should match `prop_anno`.
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
  beta = c(1, 2, 3),
  binary = FALSE,
  geno = NULL,
  include_residual = TRUE,
  indicator = FALSE,
  maf_range = c(0.001, 0.005),
  method = "none",
  n = 100,
  prop_anno = c(0.5, 0.4, 0.1),
  prop_causal = 1.0,
  random_signs = FALSE,
  random_var = 0.0, 
  snps = 100,
  weights = c(1, 1, 1)
) {

  # Generate genotypes and annotations.
  if (is.null(anno) & is.null(geno)) {
    
    # Check if annotation proportions and weights are consistent.
    if (length(prop_anno) != length(weights)) {
      stop("length(prop_anno) should equal length(weights)")
    }
    
    # Neither annotations nor genotypes provided.
    anno_geno <- GenGeno(
      n = n,
      snps = snps,
      maf_range = maf_range,
      prop_anno = prop_anno
    )
    
  } else if(is.null(anno) & !is.null(geno)) {
    
    # Only genotypes provided.
    geno <- as.matrix(geno)
    n <- nrow(geno)
    snps <- ncol(geno)
    anno <- GenAnno(snps, prop_anno = prop_anno)
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
  if (binary) {
    type <- "binary"
  } else {
    type <- "quantitative"
  }
  out <- list(
    anno = anno_geno$anno,
    covar = covar,
    geno = anno_geno$geno,
    pheno = pheno,
    type = type
  )
  return(out)
}
