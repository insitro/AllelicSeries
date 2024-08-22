# Purpose: Given simulated data from `DGP`, calculate summary statistics.
# Updated: 2024-07-31

#' Calculate Summary Statistics
#' 
#' @param anno (snps x 1) annotation vector.
#' @param covar (subjects x covars) covariate matrix.
#' @param data List of data containing the annotation vector `anno`, the 
#'   covariate data `covar`, the genotype matrix `geno`, and the phenotype
#'   vector `pheno`, as returned by \code{\link{DGP}}. Overrides the other
#'   arguments if provided.
#' @param geno (subjects x snps) genotype matrix.
#' @param pheno (subjects x 1) phenotype vector.
#' @param is_binary Is the phenotype binary? Default: FALSE.
#' 
#' @return List containing the following items: \itemize{
#' \item{anno: A SNP-length annotation vector.}
#' \item{ld: A SNP x SNP correlation (LD) matrix.}
#' \item{maf: Minor allele frequency of each variant.}
#' \item{sumstats: A SNP x 4 matrix of summary statistics.}
#' \item{type: Either "binary" or "quantitative".}
#' }
#' @export
#' @examples
#' data <- DGP()
#' sumstats <- CalcSumstats(data = data)
CalcSumstats <- function(
  anno = NULL,
  covar = NULL,
  data = NULL,
  geno = NULL,
  pheno = NULL,
  is_binary = FALSE
) {
  
  # Unpack.
  if (!is.null(data)) {
    anno <- data$anno
    covar <- data$covar
    geno <- data$geno
    pheno <- data$pheno
    is_binary <- (data$type == "binary")
  } else if (is.null(anno) || is.null(covar) || 
             is.null(geno) || is.null(pheno)) {
    
    msg <- paste0(
      "Unless 'data' as produced by 'DGP' is provided, each of ",
      "anno, covar, geno, and pheno must be provided."
    )
    stop(msg)
    
  }
  
  # Calculate MAF.
  maf <- apply(geno, 2, mean) / 2
  
  # Calculate LD matrix.
  ld <- CorCpp(geno)
  n_snp <- length(anno)
  
  # Calculate beta and SE.
  if (is_binary) {
    
    type <- "binary"
    sumstats <- lapply(seq_len(n_snp), function(i) {
      fit <- stats::glm(
        pheno ~ 0 + geno[, i] + covar,
        family = stats::binomial()
      )
      results <- summary(fit)
      beta <- as.numeric(results$coefficients[, "Estimate"][1])
      se <- as.numeric(results$coefficients[, "Std. Error"][1])
      p <- as.numeric(results$coefficients[, "Pr(>|z|)"][1])
      out <- data.frame(beta = beta, se = se, p = p)
    })
    
  } else {
    
    type <- "quantitative"
    sumstats <- lapply(seq_len(n_snp), function(i) {
      fit <- OLS(y = pheno, X = cbind(geno[, i], covar))
      beta <- fit$beta[1]
      se <- fit$se[1]
      p <- fit$pval[1]
      out <- data.frame(beta = beta, se = se, p = p)
    })
    
  }
  
  sumstats <- do.call(rbind, sumstats)
  sumstats$anno <- anno
  sumstats <- sumstats[, c("anno", "beta", "se", "p")]
  
  # Output.
  out <- list(
    anno = anno,
    ld = ld,
    maf = maf,
    sumstats = sumstats,
    type = type
  )
}

