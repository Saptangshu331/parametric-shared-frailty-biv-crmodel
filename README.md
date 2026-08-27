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
```

```html
<!--
```


## Usage

```r
library(ggflexbcr)

# data("diabetic_retinopathy")

wrk_data <- ggflexbcr::diabetic_retinopathy
attach(wrk_data)

censor_ind=rep(NA, nrow(wrk_data))

for (i in 1:nrow(wrk_data)) {

  if (status_trt[i] == 1 & status_untrt[i] == 1) {
    censor_ind[i] <- 1
  } else if (status_trt[i] == 1 & status_untrt[i] == 0) {
    censor_ind[i] <- 2
  } else if (status_trt[i] == 0 & status_untrt[i] == 1) {
    censor_ind[i] <- 3
  } else {
    censor_ind[i] <- 4
  }

}

wrk_data$censor_ind=censor_ind

age=(age-mean(age))/sd(age)

wrk_data$age=age

t1=obstime_trt
t2=obstime_untrt
censor_type_vec=censor_ind

diabetes_type=ifelse(wrk_data$type=="juvenile",1,0)

X=as.matrix(cbind(1,diabetes_type))


EM_func_qunfixed(
  c(0,0),
  c(0,0),
  .5,
  c(0,0),
  c(0,0),
  1.5,
  1.5,
  .8,
  1.2,
  t1,
  t2,
  X,
  censor_type_vec,print_opt=TRUE)


```





```html
-->
```


## Author
Saptangshu Nandi, Sandip Barui, Debanjan Mitra and Narayanaswamy Balakrishnan



## License
GPL-3 


## Notes
Users should see model assumptions before application


