
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CRTConjoint

<!-- badges: start -->

[![R-CMD-check](https://github.com/daewoongham97/CRTConjoint/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daewoongham97/CRTConjoint/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of CRTConjoint is to use the conditional randomization test
(CRT) to test for various hypothesis in conjoint experiments. In
particular, `CRT_pval` aims to test whether a factor matters in any way.
For example, does education matter in immigration preferences given
other attributes of the candidate.

## Installation

You can install CRTConjoint from
[GitHub](https://github.com/daewoongham97/CRTConjoint) with:

``` r
# install.packages("devtools")
devtools::install_github("daewoongham97/CRTConjoint")
```

or directly from CRAN with:

``` r
install.packages("CRTConjoint")
```

## Example

This is a basic example which shows you how to test whether education
matters for immigration preferences.

``` r
library(CRTConjoint)

# Immigration data
data("immigrationdata")
form = formula("Y ~ FeatEd + FeatGender + FeatCountry + FeatReason + FeatJob +
FeatExp + FeatPlans + FeatTrips + FeatLang + ppage + ppeducat + ppethm + ppgender")
left = colnames(immigrationdata)[1:9]
right = colnames(immigrationdata)[10:18]

## Not run: 
# Testing whether edcuation matters for immigration preferences
education_test = CRT_pval(formula = form, data = immigrationdata, X = "FeatEd",
 left = left, right = right, non_factor = "ppage", B = 100, analysis = 2)
education_test$p_val
```
