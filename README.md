README
================
2024-07-24

# Allelic Series

This package implements gene-level rare variant association tests
targeting allelic series: genes where increasingly deleterious mutations
have increasingly large phenotypic effects. The main COding-variant
Allelic Series Test (COAST) operates on the benign missense variants
(BMVs), deleterious missense variants (DMVs), and protein truncating
variants (PTVs) within a gene. COAST uses a set of adjustable weights
that tailor the test towards rejecting the null for genes where the
average magnitude of phenotypic effect increases monotonically from BMVs
to DMVs to PTVs. Such genes are of candidate therapeutic interest due to
the existence of a dose-response relationship between gene functionality
and phenotypic impact. See McCaw ZR, O’Dushlaine C, Somineni H, Bereket
M, Klein C, Karaletsos T, Casale FP, Koller D, Soare TW. (2023) “An
allelic-series rare-variant association test for candidate-gene
discovery”
[doi:10.1016/j.ajhg.2023.07.001](https://www.cell.com/ajhg/fulltext/S0002-9297(23)00241-0).

# Installation

``` r
install.packages("AllelicSeries")
```

``` r
library(AllelicSeries)
```

## Example data

``` r
set.seed(101)
n <- 100
data <- DGP(
  n = n,
  snps = 300,
  beta = c(1, 2, 3) / sqrt(n),
)
```

The example `data` are a list with the following components:

- `anno`: An `snps` by 1 annotation vector coded as 0 for benign
  missense variants (BMVs), 1 for deleterious missense variants (DMVs),
  and 2 for protein truncating variants (PTVs). Note that the values of
  (0, 1, 2) simply identify different categories of variants; `weights`
  other than these can be set when performing the association test.

- `covar`: An `n` by 6 covariate matrix including an intercept `int`,
  and covariates representing `age`, `sex`, and 3 genetic PCs (`pc1`,
  `pc2`, `pc3`).

- `geno`: An `n` by `snps` genotype matrix with additive coding and
  minor allele frequencies between 0.5% and 1.0%.

- `pheno`: An `n` by 1 phenotype vector.

*Note*: Scaling `beta` by `1 / sqrt(n)` makes the power invariant to the
sample size `n`.

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

The function `COAST` performs the allelic series test. The required
inputs are the annotation vector, a covariate matrix, the per-variant
genotype matrix, and the phenotype vector.

- The function assumes 3 annotation categories, coded as: `0`, `1`, `2`.
  The length of `anno` should match the number of columns of the
  genotype matrix `geno`.

- If unspecified, `covar` will default to an intercept vector (i.e. a
  vector of `1`s). If `covar` is provided, an intercept should be
  included manually, if desired.

- `weights` encodes the relative importance of BMVs, DMVs, and PTVs. The
  example weights of `c(1, 2, 3)` target a genetic architecture where
  effect sizes increase with increasing deleteriousness: BMVs have an
  effect of 1, DMVs have an effect of 2, and PTVs have an effect of 3.
  Weights of `c(1, 1, 1)` target instead a genetic architecture where
  all variant types have equivalent expected magnitudes.

``` r
show(results)
```

    ## Counts:
    ##   anno alleles variants carriers
    ## 1    0     287      163       96
    ## 2    1     162      109       82
    ## 3    2      61       28       45
    ## 
    ## 
    ## P-values:
    ##           test   type     pval
    ## 1        count burden 3.11e-26
    ## 2          ind burden 1.32e-09
    ## 3    max_count burden 3.08e-10
    ## 4      max_ind burden 5.37e-09
    ## 5    sum_count burden 1.66e-20
    ## 6      sum_ind burden 2.55e-11
    ## 7 allelic_skat   skat 2.66e-07
    ## 8         omni   omni 3.74e-25

By default, the output of `COAST` includes a data.frame of counts
showing the number of alleles, variants, and carriers in each class that
contributed to the test, and a data.frame of p-values, with the `omni`
test denoting the final, overall p-value. The counts data.frame is
accessed via:

``` r
results@Counts
```

    ##   anno alleles variants carriers
    ## 1    0     287      163       96
    ## 2    1     162      109       82
    ## 3    2      61       28       45

and the p-values data.frame via:

``` r
results@Pvals
```

    ##           test   type         pval
    ## 1        count burden 3.112702e-26
    ## 2          ind burden 1.322084e-09
    ## 3    max_count burden 3.076876e-10
    ## 4      max_ind burden 5.374363e-09
    ## 5    sum_count burden 1.661854e-20
    ## 6      sum_ind burden 2.554417e-11
    ## 7 allelic_skat   skat 2.658137e-07
    ## 8         omni   omni 3.735235e-25

To return the omnibus p-value only, specify `return_omni_only = TRUE`
when calling `COAST`.

## Robust omnibus test

In the case that all variants have comparable magnitudes, standard
SKAT-O test will have more power than an allelic series test with
weights `c(1, 2, 3)`. This is because the generative weighting scheme is
in fact `c(1, 1, 1)`. To perform a robust test that is powerful for
detecting allelic series but as powerful as standard SKAT-O when all
variants have similar magnitudes, we can incorporate standard SKAT-O in
the allelic series omnibus test. To do this, specify
`include_orig_skato_all = TRUE` when calling `COAST`. Another option is
to incorporate standard SKAT-O incorporating PTVs only:
`include_orig_skato_ptv = TRUE`.

``` r
results <- COAST(
  anno = data$anno,
  covar = data$covar,
  geno = data$geno,
  pheno = data$pheno,
  include_orig_skato_all = TRUE,
  include_orig_skato_ptv = TRUE,
  weights = c(1, 2, 3)
)
show(results)
```

    ## Counts:
    ##   anno alleles variants carriers
    ## 1    0     287      163       96
    ## 2    1     162      109       82
    ## 3    2      61       28       45
    ## 
    ## 
    ## P-values:
    ##             test   type     pval
    ## 1          count burden 3.11e-26
    ## 2            ind burden 1.32e-09
    ## 3      max_count burden 3.08e-10
    ## 4        max_ind burden 5.37e-09
    ## 5      sum_count burden 1.66e-20
    ## 6        sum_ind burden 2.55e-11
    ## 7   allelic_skat   skat 2.66e-07
    ## 8  orig_skat_all   skat 1.55e-05
    ## 9  orig_skat_ptv   skat 6.63e-08
    ## 10          omni   omni 3.74e-25

## Loading genotypes

The [genio](https://CRAN.R-project.org/package=genio) and
[rbgen](https://enkre.net/cgi-bin/code/bgen/wiki?name=rbgen) packages
may be used to load PLINK and BGEN genotypes in R, respectively.
Moreover, [PLINK](https://www.cog-genomics.org/plink/2.0/) enables
conversion between these file types.
