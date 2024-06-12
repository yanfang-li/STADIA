
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STADIA <img src="man/figures/logo.png" align="right" width="160" />

<!--[![R-CMD-check](https://github.com/yanfang-li/stadia/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yanfang-li/stadia/actions/workflows/R-CMD-check.yaml)
<!-- [![Codecov test coverage](https://codecov.io/gh/yanfang-li/stadia/branch/master/graph/badge.svg)](https://app.codecov.io/gh/yanfang-li/stadia?branch=master) -->

[![CRAN
Version](https://www.r-pkg.org/badges/version/stadia)](https://cran.r-project.org/package=stadia)

## Overview

STADIA is designed to integrate multiple spatial transcriptomics slices,
which simultaneously achieve

- dimension reduction
- batch effects correction
- spatial domains identification

Overall, STADIA is a Bayesian hierarchical hidden Markov random field
model that first uses Bayesian factor regression and location-and-scale
adjustment to learn a batch-corrected low-dimensional representation of
the gene expression profiles, and then spatially clusters the embedding
using a Gaussian mixture model with a Potts spatial prior to promote
local consistency.

## Installation

<!-- Install the released version of stadia R package from CRAN:
&#10;
```r
# Install released version from CRAN
install.packages("stadia")
```
&#10;or install the development version from GitHub with:-->

It’s recommended to create a separate `conda` environment for running
STADIA:

``` r
conda create -n env_stadia -c conda-forge r-base r-essentials
conda activate env_stadia
```

First, make sure all dependencies are available:

``` r
dependencies <- c('Seurat','mclust','irlba','mombf','progress','Rcpp','BiocSingular','BiocParallel','BiocNeighbors','RcppArmadillo','RcppDist')

ap <- available.packages()
install.packages(dependencies[!(dependencies %in% row.names(installed.packages())) & (dependencies %in% row.names(ap))])

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(dependencies[!(dependencies %in% row.names(installed.packages())) & !(dependencies %in% row.names(ap))])
```

Install the stadia package from GitHub:

``` r
# Install released version from github
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yanfang-li/stadia")
```

or from the source in the terminal:

``` r
R CMD INSTALL stadia_1.0.0.tar.gz
```

The stadia package has been built and tested on the following operating
system:

- Windows: 10, version 22H2 \[Rtools required\]
- Linux: Ubuntu 20.04.4 LTS
- MacOS: Monterey 12.7.1 (Apple M1), Monterey 12.0.1 (Intel Core i7)

## Setup to install the stadia package on MacOS.

Since the stadia package uses C++ and openmp, some additional
configuration is required to install it from source on MacOS. For
details, please refer to the following steps:

1.  Install **Command Line Tools for Xcode** from
    <https://developer.apple.com/download/all>.
2.  Download **homebrew** from
    <https://github.com/Homebrew/brew/releases/tag/4.1.24> and install
    it.
3.  Install gcc using **brew install gcc**.
4.  Install LLVM using **brew install llvm**.
5.  Install libomp package using **brew install libomp**.
6.  Make file **~/.R/Makevars**: (all paths must to be replaced with the
    paths on your own computer)

``` bash
# MacOS Xcode header location
# (do "xcrun -show-sdk-path" in terminal to get path)
XC_LOC = /Applications/Xcode.13.4.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
# (do "echo $(brew --prefix **)" in terminal to get path)
LLVM_LOC = /opt/homebrew/opt/llvm
GCC_LOC = /opt/homebrew/Cellar/gcc/12.2.0
GETTEXT_LOC = /opt/homebrew/opt/gettext
OMP_LOC= /opt/homebrew/opt/libomp
OPENSSL_LOC=/opt/homebrew/opt/openssl@1.1

CC=$(LLVM_LOC)/bin/clang
CXX=$(LLVM_LOC)/bin/clang++
CXX11=$(LLVM_LOC)/bin/clang++
CXX14=$(LLVM_LOC)/bin/clang++
CXX17=$(LLVM_LOC)/bin/clang++
CXX1X=$(LLVM_LOC)/bin/clang++

CFLAGS=-g -O2 -Wall -pedantic -std=gnu99 -mtune=native -pipe
CXXFLAGS=-g -O2 -Wall -pedantic -std=c++11 -mtune=native -pipe
LDFLAGS=-L"$(LLVM_LOC)/lib" -L"$(GETTEXT_LOC)/lib" -Wl,-rpath,$(LLVM_LOC)/lib --sysroot="$(XC_LOC)" -lomp -L"$(OPENSSL_LOC)/lib"
CPPFLAGS=-I"$(GETTEXT_LOC)/include" -I"$(LLVM_LOC)/include" -I"$(OMP_LOC)/include" -isysroot "$(XC_LOC)" -Xclang -fopenmp -I"$(OPENSSL_LOC)/include"

FC=$(GCC_LOC)/bin/gfortran
F77=$(GCC_LOC)/bin/gfortran
FLIBS=-L$(GCC_LOC)/lib/gcc/11/ -lm
```

## Steps to use

The main steps to run the STADIA algorithm using the *stadia* package
are (it is recommended to run in the server background for large data
sets)

``` r
## load packages
library(stadia)
## set hyperparameters
hyper <- HyperParameters(obj, dim = d, eta = eta)
## run model
out <- stadia(obj, hyper, dim = d, n_cluster = K, 
             platform = "visium")
```

where

- obj: a list of Seurat objects containing information of ST slices
- dim: the dimension of latent factors
- eta: the spatial smoothness parameter used in Potts model
- n_cluster: the number of spatial domains
- platform: specifics the platform

## Demonstration

Files in the
[Applications](https://yanfang-li.github.io/stadia/articles/stadia.html)
demonstrate how to use the stadia package to run the STADIA algorithm in
the `Run` section.

## Paper Citation

Yanfang Li and Shihua Zhang (2023+). “Statistical batch-aware embedded
integration, dimension reduction and alignment for spatial
transcriptomics”.
