# Allelic Series

Implementation of gene-level rare coding variant association tests
targeting allelic series: cases where increasingly deleterious mutations
have increasingly large phenotypic effects.

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
  and 2 for protein truncating variants (PTVs).

- `covar`: An `n` by 6 covariate matrix including an intercept `int`,
  and covariates representing `age`, `sex`, and 3 genetic PCs (`pc1` …
  `pc3`).

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
  BMVs have no effect and PTVs have twice the effect of DMVs. Weights of
  `c(1, 1, 1)` target instead a genetic architecture where all variant
  types have equivalent expected magnitudes.

``` r
show(results)
```

    ##        p_count          p_ind    p_max_count      p_max_ind    p_sum_count 
    ##   7.707024e-29   6.745269e-06   4.299938e-13   6.228669e-07   6.756953e-18 
    ##      p_sum_ind p_allelic_skat         p_omni 
    ##   1.894500e-06   8.507977e-08   9.248429e-28

By default, the output of `COAST` is a vector of p-values, corresponding
to the different components of the allelic series test and the overall
omnibus test (`p_omni`). To return the omnibus p-value only, specify
`return_omni_only = TRUE` when calling `COAST`.

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

    ##         p_count           p_ind     p_max_count       p_max_ind     p_sum_count 
    ##    7.707024e-29    6.745269e-06    4.299938e-13    6.228669e-07    6.756953e-18 
    ##       p_sum_ind  p_allelic_skat p_orig_skat_all p_orig_skat_ptv          p_omni 
    ##    1.894500e-06    8.507977e-08    1.250426e-05    2.237988e-13    9.248429e-28
