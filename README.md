# Parametric Shared Frailty Bivariate Cure Rate Model

## Overview

This R package implements a parametric bivariate cure rate model with shared frailty using Rcpp for efficient computation.

It is designed for modeling:
- Dependent survival outcomes
- Presence of cure combinations
- Effect of covariates on lifetimes and cure combinations

---


## Model Description

Refer to paper


## Package Structure


R/ — R interface functions, mainly EM algorithm

src/ — Contains C++ implementation using Rcpp and related libraries for faster computation. The functions support the EM algorithm.

DESCRIPTION — Package metadata

NAMESPACE — Exported functions and imports from other packages.

DESCRIPTION
Provides package metadata such as name, version, dependencies, author information, and licensing.
NAMESPACE
Defines exported functions and manages imports from other packages.


## Installation

You can install the package **ggflexbcr** from GitHub:

```r
# install.packages("remotes")  # if not already installed
remotes::install_github("Saptangshu331/parametric-shared-frailty-biv-crmodel")

## Author
Saptangshu Nandi, Sandip Barui, Debanjan Mitra and Narayanaswamy Balakrishnan



## License
GPL-3 


## Notes
Users should see model assumptions before application


