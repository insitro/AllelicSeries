# Purpose: Output class.
# Updated: 2024-11-20

#' Data.frame or Null Class
#' 
#' @name DfOrNULL-class
#' @rdname DfOrNULL-class
#' @importClassesFrom methods data.frame
#' @exportClass DfOrNULL
setClassUnion("DfOrNULL", members = c("data.frame", "NULL"))


#' Allelic Series Output Class
#' 
#' @slot Betas Effect sizes and standard errors.
#' @slot Counts Allele, variant, and carrier counts.
#' @slot Pvals Result p-values.
#' @name COAST-class
#' @rdname COAST-class
#' @exportClass COAST
setClass(
  Class = "COAST",
  representation = representation(
    Betas = "DfOrNULL",
    Counts = "DfOrNULL",
    Pvals = "data.frame"
  )
)


#' Print Method for COAST Object.
#'
#' Print method for objects of class \code{COAST}.
#'
#' @param x An object of class \code{COAST}.
#' @param ... Unused.
#' @export
print.COAST <- function(x, ...) {

  # Betas
  if (!is.null(x@Betas)) {
    cat("Effect Sizes:\n")
    betas <- x@Betas
    betas$beta <- formatC(betas$beta, format = "f", digits = 2)
    betas$se <- formatC(betas$se, format = "f", digits = 3)
    show(betas)
    cat("\n\n")
  }
  
  # Counts.
  if (!is.null(x@Counts)) {
    cat("Counts:\n")
    show(x@Counts)
    cat("\n\n")
  }
  
  # P-values.
  cat("P-values:\n")
  pvals <- x@Pvals
  pvals$pval <- formatC(pvals$pval, format = "e", digits = 2)
  show(pvals)
  cat("\n\n")
}


#' Show Method for COAST Object
#'
#' @param object An object of class \code{COAST}.
#' @rdname COAST-method
#' @importFrom methods show
setMethod(
  f = "show",
  signature = c(object = "COAST"),
  definition = function(object) {
    print.COAST(x = object)
  }
)

