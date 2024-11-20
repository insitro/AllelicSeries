README
================

Updated: 2024-11-20

# Allelic Series

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/AllelicSeries)](https://cran.r-project.org/package=AllelicSeries)
[![](https://cranlogs.r-pkg.org/badges/grand-total/AllelicSeries)](https://CRAN.R-project.org/package=AllelicSeries)

This package implements gene-level rare variant association tests
targeting allelic series: genes where increasingly deleterious mutations
have increasingly large phenotypic effects. The main COding-variant
Allelic Series Test (COAST) was designed to operate on the benign
missense variants (BMVs), damaging missense variants (DMVs), and protein
truncating variants (PTVs) within a gene. COAST uses a set of adjustable
weights that tailor the test towards rejecting the null for genes where
the mean phenotypic impact increases monotonically from BMVs to DMVs to
PTVs. Such genes are of candidate therapeutic interest due to the
existence of a dose-response relationship between gene functionality and
phenotypic effect. COAST-SS extends COAST to accept summary statistics
as input. Both methods have been generalized to allow for arbitrary
number of discrete annotation categories, for example the 4 category
scale: benign missense, damaging missense, low-confidence loss of
function, high-confidence loss of function. For additional details, see:

- McCaw ZR, O’Dushlaine C, Somineni H, Bereket M, Klein C, Karaletsos T,
  Casale FP, Koller D, Soare TW. (2023) “An allelic-series rare-variant
  association test for candidate-gene discovery”
  [doi:10.1016/j.ajhg.2023.07.001](https://www.cell.com/ajhg/fulltext/S0002-9297(23)00241-0).

- McCaw ZR, Gao J, Dey R, Research Team insitro, Fox E, O’Dushlaine C,
  Soare TW. (2024) “Extending the Coding-variant Allelic Series Test to
  Summary Statistics”
  [doi:10.1101/2024.10.31.621375](https://www.biorxiv.org/content/10.1101/2024.10.31.621375v1).

# Installation

``` r
# Install from CRAN.
install.packages("AllelicSeries")

# Install from GitHub.
devtools::install_github(repo = "insitro/AllelicSeries")
```

``` r
library(AllelicSeries)
```

To view vignettes:

``` r
browseVignettes("AllelicSeries")
```

## Example data

Here, data for `n = 1000` subjects at `snps = 300` variants with minor
allele frequencies between 0.1% and 0.5% (`maf_range`) are simulated. By
default, variants are annotated into 1 of 3 categories, representing
BMVs, DMVs, and PTVs. `beta` specifies the mean per-variant effect sizes
in each annotation category. Setting `beta` inversely proportional to
`sqrt(n)` makes the power independent of the sample size. For additional
details, see the Data Generation vignette.

``` r
set.seed(101)
n <- 1000
data <- DGP(
  n = n,
  snps = 300,
  maf_range = c(0.001, 0.005),
  beta = c(1, 2, 3) / sqrt(n)
)
```

The example `data` are a list with the following components:

- `anno`: An `snps` by 1 annotation vector coded as 1 for BMVs, 2 for
  DMVs, and 3 for PTVs. Note that the values of (1, 2, 3) simply
  identify different categories of variants; `weights` other than these
  can be set when performing the association test.

- `covar`: An `n` by 6 covariate matrix where `n` is the number of
  subjects and the columns correspond to: an intercept `int`, `age`,
  `sex`, and 3 genetic PCs (`pc1`, `pc2`, `pc3`).

- `geno`: An `n` by `snps` genotype matrix with additive coding.

- `pheno`: An `n` by 1 phenotype vector.

## COding-variant Allelic Series Test

``` r
results <- COAST(
  anno = data$anno,
  covar = data$covar,
  geno = data$geno,
  pheno = data$pheno,
  weights = c(1, 2, 3)
)
```

The function `COAST` performs the coding-variant allelic series test.
The required inputs are the annotation vector, a covariate matrix, the
per-variant genotype matrix, and the phenotype vector.

- The function allows any number of annotation categories, coded as
  consecutive integers starting at 1. The example data have categories
  of `c(1, 2, 3)`. The length of `anno` should match the number of
  columns of the genotype matrix `geno`. Note: a previous version of the
  package used zero-indexed annotation categories. The main functions
  will still run with zero-indexed annotations, but one-based indexing
  has been adopted as it is the default in `R`.

- If unspecified, `covar` will default to an intercept vector (i.e. a
  vector of `1`s). If `covar` is provided, an intercept should be
  included manually as necessary.

- `weights` encodes the relative importance of the annotation
  categories. The example weights of `c(1, 2, 3)` target a genetic
  architecture where effect sizes increase with increasing
  deleteriousness: BMVs have a relative effect of 1, DMVs have a
  relative effect of 2, and PTVs have a relative effect of 3. Weights of
  `c(1, 1, 1)` target instead a genetic architecture where all variant
  types have equivalent expected magnitudes.

``` r
show(results)
```

    ## Effect Sizes:
    ##         test beta    se
    ## 1       base 0.03 0.027
    ## 2       base 0.07 0.028
    ## 3       base 0.07 0.062
    ## 4        ind 0.01 0.048
    ## 5        ind 0.16 0.048
    ## 6        ind 0.05 0.068
    ## 7  max_count 0.03 0.015
    ## 8    max_ind 0.07 0.025
    ## 9  sum_count 0.03 0.010
    ## 10   sum_ind 0.04 0.016
    ## 
    ## 
    ## Counts:
    ##   anno alleles variants carriers
    ## 1    1     834      148      574
    ## 2    2     730      126      516
    ## 3    3     154       26      144
    ## 
    ## 
    ## P-values:
    ##           test   type     pval
    ## 1     baseline burden 2.64e-02
    ## 2          ind burden 5.84e-03
    ## 3    max_count burden 2.34e-02
    ## 4      max_ind burden 8.95e-03
    ## 5    sum_count burden 2.62e-03
    ## 6      sum_ind burden 4.39e-03
    ## 7 allelic_skat   skat 3.38e-03
    ## 8         omni   omni 4.36e-03

By default, the output of `COAST` includes a data.frame of estimated
effect sizes from the burden tests, a data frame of counts showing the
number of alleles, variants, and carriers in each class that contributed
to the test, and a data.frame of p-values, with the `omni` test denoting
the final, overall p-value. The effect sizes data.frame is accessed via:

``` r
results@Betas
```

    ##         test        beta         se
    ## 1       base 0.029482960 0.02667753
    ## 2       base 0.070413431 0.02760127
    ## 3       base 0.073294520 0.06157383
    ## 4        ind 0.007297359 0.04823273
    ## 5        ind 0.163297146 0.04755333
    ## 6        ind 0.054511276 0.06780962
    ## 7  max_count 0.034206336 0.01508838
    ## 8    max_ind 0.066202953 0.02532496
    ## 9  sum_count 0.031465473 0.01045514
    ## 10   sum_ind 0.044181239 0.01551063

the counts data.frame is accessed via:

``` r
results@Counts
```

    ##   anno alleles variants carriers
    ## 1    1     834      148      574
    ## 2    2     730      126      516
    ## 3    3     154       26      144

and the p-values data.frame via:

``` r
results@Pvals
```

    ##           test   type        pval
    ## 1     baseline burden 0.026375086
    ## 2          ind burden 0.005839884
    ## 3    max_count burden 0.023386302
    ## 4      max_ind burden 0.008945286
    ## 5    sum_count burden 0.002616191
    ## 6      sum_ind burden 0.004393287
    ## 7 allelic_skat   skat 0.003376729
    ## 8         omni   omni 0.004363169

To return the omnibus p-value only, specify `return_omni_only = TRUE`
when calling `COAST`. For additional details, see the COAST vignette.

## Different numbers of annotation categories

The following provides a minimal example of generating data and running
COAST with a different number of annotation categories, in this case 4,
which might represent benign missense variants, damaging missense
variants, low-confidence loss of function, and high-confidence loss of
function. The main difference when analyzing data with a different
number of annotation categories is that the `weights` argument should be
provided to `COAST`, and should have length equal to the number of
annotation categories.

``` r
withr::local_seed(102)

# Generate data.
n <- 1e2
data4 <- DGP(
  n = n,
  snps = 400,
  beta = c(1, 2, 3, 4) / sqrt(n),
  prop_anno = c(0.4, 0.3, 0.2, 0.1),
  weights = c(1, 1, 1, 1)
)

# Run COAST-SS.
results <- COAST(
  anno = data4$anno,
  covar = data4$covar,
  geno = data4$geno,
  pheno = data4$pheno,
  weights = c(1, 2, 3, 4)
)
show(results)
```

    ## Effect Sizes:
    ##         test  beta    se
    ## 1       base  0.01 0.065
    ## 2       base  0.23 0.062
    ## 3       base  0.19 0.096
    ## 4       base  0.61 0.116
    ## 5        ind -0.10 0.253
    ## 6        ind  0.53 0.205
    ## 7        ind  0.31 0.169
    ## 8        ind  0.96 0.186
    ## 9  max_count  0.18 0.035
    ## 10   max_ind  0.47 0.090
    ## 11 sum_count  0.11 0.018
    ## 12   sum_ind  0.19 0.035
    ## 
    ## 
    ## Counts:
    ##   anno alleles variants carriers
    ## 1    1     194      169       86
    ## 2    2     159      131       78
    ## 3    3      71       61       51
    ## 4    4      41       39       31
    ## 
    ## 
    ## P-values:
    ##           test   type     pval
    ## 1     baseline burden 9.67e-10
    ## 2          ind burden 5.86e-07
    ## 3    max_count burden 4.26e-07
    ## 4      max_ind burden 1.65e-07
    ## 5    sum_count burden 5.37e-10
    ## 6      sum_ind burden 8.24e-08
    ## 7 allelic_skat   skat 3.27e-06
    ## 8         omni   omni 4.11e-09

## COAST from Summary Statistics

### Summary statistics calculation

The function `CalcSumstats` calculates summary statistics starting
either from:

- The direct output of `DGP`, or
- An annotation vector `anno`, genotype matrix `geno`, and phenotype
  vector `pheno`, formatted as provided by `DGP`.
- Providing covariates `covar` is optional. If omitted, an
  intercept-only design matrix is adopted by default. If supplied, the
  covariates should include an intercept as necessary.

``` r
sumstats <- CalcSumstats(
  anno = data$anno,
  covar = data$covar,
  geno = data$geno,
  pheno = data$pheno
)
```

The output `sumstats` is a list containing:

- `ld`, a (snps x snps) LD (genotype correlation) matrix.
- `sumstats`, a (snps x 5) data.frame including the annotations `anno`,
  minor allele frequencies `maf`, effect sizes `beta`, standard errors
  `se`, and p-values `p`.

``` r
head(sumstats$sumstats)
```

    ##   anno   maf        beta        se         p
    ## 1    1 4e-03 -0.06655739 0.2688400 0.8044652
    ## 2    2 4e-03 -0.01192392 0.2685362 0.9645829
    ## 3    1 1e-03 -0.44322853 0.5355026 0.4078478
    ## 4    1 5e-04  0.41017824 0.7583407 0.5885840
    ## 5    1 5e-03 -0.01965152 0.2417292 0.9352069
    ## 6    1 3e-03 -0.07339959 0.3096588 0.8126306

### Running COAST from summary statistics

`COASTSS` is the main function for running the coding-variant allelic
series test from summary statistics. The necessary inputs are the
annotation vector `anno`, the effect size vector `beta`, and the
standard error vector `se`. The test will run *without* the LD matrix
`ld` matrix. However, to do so, it assumes the variants are in linkage
equilibrium (i.e. `ld` is the identity matrix). This approximation is
reasonable when the LD is minimal, as is expected among rare variants,
however it may break down if variants of sufficient minor allele count
are included in the analysis. If available, providing the in-sample LD
matrix is always recommended. If provided, the `maf` vector is used to
up-weight rare variants during the allelic SKAT test.

``` r
results <- COASTSS(
  anno = sumstats$sumstats$anno,
  beta = sumstats$sumstats$beta,
  se = sumstats$sumstats$se,
  maf = sumstats$sumstats$maf,
  ld = sumstats$ld,
  weights = c(1, 2, 3)
)
show(results)
```

    ## P-values:
    ##           test   type     pval
    ## 1     baseline burden 2.75e-02
    ## 2    sum_count burden 2.63e-03
    ## 3 allelic_skat   skat 3.23e-03
    ## 4         omni   omni 3.86e-03

In comparing the outputs of the summary statistics based test to those
of the individual level data test, several differences are noteworthy:

- Not all components of `COAST` could be included in `COASTSS`. In
  particular, the `max` tests cannot be obtained starting from standard
  summary statistics. In addition, by convention, summary statistics are
  generated from count rather than indicator genotypes. If available,
  `COASTSS` can be applied to summary statistics generated from
  indicator genotypes.

- Several approximations are required in order to perform the
  coding-variant allelic series test with summary statistics. As such,
  the p-values obtained from `COASTSS` and `COAST` will not be
  identical, even when starting from the same data. Nonetheless, the
  operating characteristics of `COASTSS` (and the original `COAST`) have
  been validated through extensive simulation studies.

- `COASTSS` can also be run with a different number of annotation
  categories. To do so, specify a `weight` vector of the appropriate
  length.

For additional details, see the COAST-SS vignette.

# Appendix

## Loading genotypes

The [genio](https://CRAN.R-project.org/package=genio) and
[rbgen](https://enkre.net/cgi-bin/code/bgen/wiki?name=rbgen) packages
may be used to load PLINK and BGEN genotypes in R, respectively.
Moreover, [PLINK](https://www.cog-genomics.org/plink/2.0/) enables
conversion between these file types.
