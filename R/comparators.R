# Purpose: Comparators for allelic series test.
# Updated: 2023-10-04


#' Comparator Test
#'
#' Runs burden, SKAT, and SKAT-O, using default settings.
#'
#' @param covar (n x p) covariate matrix.
#' @param geno (n x snps) genotype matrix.
#' @param pheno (n x 1) phenotype vector.
#' @param apply_int Apply rank-based inverse normal transform to the phenotype?
#'   Default: TRUE. Ignored if phenotype is binary.
#' @param is_pheno_binary Is the phenotype binary? Default: FALSE.
#' @return Numeric vector of p-values.
#' @examples 
#' # Generate data.
#' data <- DGP(n = 1e3, snps = 1e2)
#' 
#' # Run the comparators.
#' results <- Comparator(
#'   geno = data$geno,
#'   pheno = data$pheno,
#'   covar = data$covar
#' )
#' @export
Comparator <- function(
    covar,
    geno,
    pheno,
    apply_int = TRUE,
    is_pheno_binary = FALSE
) {

  # Rank-normal phenotype.
  if (!is_pheno_binary & apply_int) {
    y <- RNOmni::RankNorm(pheno)
  } else {
    y <- pheno
  }

  # SKAT null model.
  if (is_pheno_binary) {
    out_type <- "D"
  } else {
    out_type <- "C"
  }
  null <- SKAT::SKAT_Null_Model(y ~ covar, out_type = out_type)

  # Alternate allele frequencies.
  aaf <- apply(geno, 2, mean) / 2
  v <- aaf * (1 - aaf)
  skat_weights <- sqrt(1 / v)

  # SKAT test.
  burden_test <- SKAT::SKAT(geno, null, r.corr = 1.0)
  p_burden <- burden_test$p.value

  skat_test <- SKAT::SKAT(
    geno,
    null,
    r.corr = 0.0,
    weights = skat_weights
  )
  p_skat <- skat_test$p.value

  skato_test <- SKAT::SKAT(
    geno,
    null,
    method = "SKATO",
    weights = skat_weights
  )
  p_skato <- skato_test$p.value

  # Output.
  out <- c(
    p_burden = p_burden,
    p_orig_skat = p_skat,
    p_orig_skato = p_skato
  )

  return(out)
}
