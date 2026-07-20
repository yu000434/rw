
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `R/rw`: Robins-Wang variance estimation for multiple imputation

<!-- badges: start -->

[![R-CMD-check](https://github.com/LucyMcGowan/rw/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucyMcGowan/rw/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/LucyMcGowan/rw/graph/badge.svg)](https://app.codecov.io/gh/LucyMcGowan/rw)
<!-- badges: end -->

The purpose of the `rw` package is to compute [Robins-Wang variance
estimates](https://doi.org/10.1093/biomet/87.1.113) for multiply imputed
data analyses. The package implements RW-MICE for parametric
chained-equation imputation using `method = "norm"` or
`method = "logreg"`, and RW-PMM for predictive mean matching using
`method = "pmmrw"`.

## Installation

This package requires the development version of mice with the `tasks`
argument:

``` r
remotes::install_github("amices/mice@dev")
```

Then you can install the development version of `rw` like so:

``` r
remotes::install_github("LucyMcGowan/rw")
```

## RW-MICE with parametric imputation

The first step is to impute your data using `mice`. When using the
`mice` function, you must set the parameter `tasks = "train"` in order
to pass the necessary parts to our subsequent functions. In this
example, we use the `nhanes` data from the mice package and impute the
incomplete variables `bmi`, `hyp`, and `chl` using normal working
models. The variable `age` has no missing values.

``` r
library(rw)
library(mice)
set.seed(1)
nhanes_scaled <- as.data.frame(scale(nhanes))

meth_norm <- make.method(nhanes_scaled)
meth_norm[] <- ""
meth_norm[c("bmi", "hyp", "chl")] <- "norm"

imp_norm <- mice(nhanes_scaled, 
                 method = meth_norm,
                 m = 5, 
                 tasks = "train",
                 print = FALSE)
```

Now, suppose we want to fit a model predicting `bmi` from `age` and
`hyp`. The `with_rw` function allows this, where the first argument is
the imputation object from the `mice` function above and the second is
the model expression. Note that currently outcome models must be
Gaussian or Binomial.

``` r
fit_rw_norm <- with_rw(imp_norm, lm(bmi ~ age + hyp))
```

Finally, we can pool the results from the fit object and calculate the
Robins-Wang variance using the `pool_rw` function

``` r
pool_rw(fit_rw_norm)
#> 
#> ── Robins-Wang Pooled Results ──────────────────────────────────────────────────
#> Number of imputations: 5
#> Sample size: 25
#> 
#>                    term estimate std.error statistic p.value conf.low conf.high
#> (Intercept) (Intercept) -0.02128    0.7101  -0.02997  0.9764   -1.494     1.451
#> age                 age -0.68330    1.5840  -0.43137  0.6704   -3.968     2.602
#> hyp                 hyp  0.15922    2.1031   0.07570  0.9403   -4.202     4.521
```

Let’s compare this result to using Rubin’s rules.

``` r
fit_rr_norm <- with(imp_norm, lm(bmi ~ age + hyp))
pool(fit_rr_norm) |>
  summary()
#>          term    estimate std.error   statistic        df    p.value
#> 1 (Intercept) -0.02127991 0.2201983 -0.09663976 16.243158 0.92419483
#> 2         age -0.68329820 0.2886778 -2.36699283 12.046683 0.03552134
#> 3         hyp  0.15921513 0.3108677  0.51216365  5.073478 0.63004620
```

## RW-PMM with predictive mean matching

RW-PMM is the Robins-Wang variance estimator for analyses after
predictive mean matching. The `pmmrw` method produces the same
imputations as standard `mice` PMM while storing the additional
information required by `pool_rw()`.

Here we impute `bmi` using `age`, then analyze `bmi ~ age`.

``` r
nhanes_pmm <- nhanes_scaled[c("age", "bmi")]

meth_pmm <- make.method(nhanes_pmm)
meth_pmm[] <- ""
meth_pmm["bmi"] <- "pmmrw"

set.seed(1)
imp_pmm <- mice(nhanes_pmm,
                method = meth_pmm,
                m = 5,
                tasks = "train",
                print = FALSE)
```

We can verify that `pmmrw` does not change the imputed values.

``` r
meth_standard_pmm <- meth_pmm
meth_standard_pmm["bmi"] <- "pmm"

set.seed(1)
imp_standard_pmm <- mice(nhanes_pmm,
                         method = meth_standard_pmm,
                         m = 5,
                         tasks = "train",
                         print = FALSE)

identical(imp_standard_pmm$imp$bmi, imp_pmm$imp$bmi)
#> [1] TRUE
```

``` r
fit_rw_pmm <- with_rw(imp_pmm, lm(bmi ~ age))
pool_rw(fit_rw_pmm)
#> 
#> ── Robins-Wang Pooled Results ──────────────────────────────────────────────────
#> Number of imputations: 5
#> Sample size: 25
#> 
#>                    term estimate std.error statistic p.value conf.low conf.high
#> (Intercept) (Intercept)  0.06773    0.3036    0.2231  0.8254  -0.5602    0.6957
#> age                 age -0.43378    0.2968   -1.4617  0.1573  -1.0477    0.1801
```

Compare this with Rubin’s rules for the same PMM imputations.

``` r
fit_rr_pmm <- with(imp_pmm, lm(bmi ~ age))
pool(fit_rr_pmm) |>
  summary()
#>          term    estimate std.error  statistic        df    p.value
#> 1 (Intercept)  0.06773122 0.2362112  0.2867401  8.375472 0.78127746
#> 2         age -0.43377882 0.2105222 -2.0604897 13.913963 0.05856272
```

## Giganti & Shepherd Example

Below is an example to replicate the [Giganti & Shepherd
(2020)](https://doi.org/10.1093/aje/kwaa153) results. This example uses
the parametric `norm` and `logreg` imputation models and demonstrates
the original Robins-Wang workflow; it is not a PMM correlation-adjusted
example.

``` r
set.seed(1)
meth <- make.method(giganti_data)
meth[] <- ""
meth[c("A", "D")] <- c("norm", "logreg")

pred <- make.predictorMatrix(giganti_data)
pred[,] <- 0
pred["A", c("X1", "X2", "A.star", "D.star")] <- 1
pred["D", c("X1", "X2", "A.star", "D.star", "A")] <- 1

imp <- mice(giganti_data, m = 10, method = meth, predictorMatrix = pred,
            print = FALSE, tasks = "train")

fit_rw <- with_rw(imp, glm(D ~ A, subset = A > 2, family = binomial()))

pooled <- pool_rw(fit_rw)
pooled
#> 
#> ── Robins-Wang Pooled Results ──────────────────────────────────────────────────
#> Number of imputations: 10
#> Sample size: 4000
#> 
#>                    term estimate std.error statistic   p.value conf.low
#> (Intercept) (Intercept)  -3.4864   0.23876    -14.60 2.721e-48  -3.9543
#> A                     A   0.8764   0.08217     10.67 1.469e-26   0.7154
#>             conf.high
#> (Intercept)    -3.018
#> A               1.038
```
