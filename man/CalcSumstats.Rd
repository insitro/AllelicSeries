% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_sumstats.R
\name{CalcSumstats}
\alias{CalcSumstats}
\title{Calculate Summary Statistics}
\usage{
CalcSumstats(
  anno = NULL,
  covar = NULL,
  data = NULL,
  geno = NULL,
  pheno = NULL,
  is_binary = FALSE
)
}
\arguments{
\item{anno}{(snps x 1) annotation vector.}

\item{covar}{(subjects x covars) covariate matrix.}

\item{data}{List of data containing the annotation vector \code{anno}, the
covariate data \code{covar}, the genotype matrix \code{geno}, and the phenotype
vector \code{pheno}, as returned by \code{\link{DGP}}. Overrides the other
arguments if provided.}

\item{geno}{(subjects x snps) genotype matrix.}

\item{pheno}{(subjects x 1) phenotype vector.}

\item{is_binary}{Is the phenotype binary? Default: FALSE.}
}
\value{
List containing the following items: \itemize{
\item{anno: A SNP-length annotation vector.}
\item{ld: A SNP x SNP correlation (LD) matrix.}
\item{maf: Minor allele frequency of each variant.}
\item{sumstats: A SNP x 4 matrix of summary statistics.}
\item{type: Either "binary" or "quantitative".}
}
}
\description{
Calculate Summary Statistics
}
\examples{
data <- DGP()
sumstats <- CalcSumstats(data = data)
}
