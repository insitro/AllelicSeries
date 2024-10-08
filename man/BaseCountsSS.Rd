% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{BaseCountsSS}
\alias{BaseCountsSS}
\title{Baseline Counts Test from Sumstats}
\usage{
BaseCountsSS(beta, ld, se)
}
\arguments{
\item{beta}{(snps x 1) vector of effect sizes for
the coding genetic variants within a gene.}

\item{ld}{(snps x snps) matrix of correlations among the genetic variants.}

\item{se}{(snps x 1) vector of standard errors for the effect sizes.}
}
\value{
Numeric p-value.
}
\description{
Baseline Counts Test from Sumstats
}
