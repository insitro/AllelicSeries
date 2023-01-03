#' @useDynLib AllelicSeries, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' Allelic Series Test
#' 
#' Implementation of gene-level rare variant association tests targeting allelic
#' series: genes where increasingly deleterious mutations have increasingly
#' large phenotypic effects. The COding-variant Allelic Series Test (COAST)
#' operates on the benign missense variants (BMVs), deleterious missense
#' variants (DMVs), and protein truncating variants (PTVs) within a gene. COAST
#' uses a set of adjustable weights that tailor the test towards rejecting the
#' null hypothesis for genes where the average magnitude of effect increases
#' monotonically from BMVs to DMVs to PTVs. See McCaw ZR, Somineni H, Bereket M,
#' Klein C, Karaletsos T, Casale FP, Koller D, Soare TW. "An allelic series rare
#' variant association test for candidate gene discovery"
#' <doi:10.1101/2022.12.23.521658>.
#' 
#' @author Zachary R. McCaw, Christoph Klein
#' @docType package
#' @name AllelicSeries-help
NULL