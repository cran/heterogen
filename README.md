# heterogen: An R package with spatial functions for heterogeneity and climate variability

P. Joser Atauchi & A. Townsend Peterson

## Installation

`heterogen` is available from CRAN, so you can use `install.packages("heteogen")` to get the current *released version*. Current version was tested on Windows 11 and OSX 15+.

The easiest way to use the *development version* on Windows or MacOS, is to install it from github:

### From source-code

To install from source-code, first install the [Rcpp](https://cran.r-project.org/package=Rcpp), [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo), and [RcppEigen](https://cran.r-project.org/package=RcppEigen) package that `heterogen` depends on:

```         
install.packages(c("Rcpp","RcppArmadillo","RcppEigen"))
```

And then continue based on the OS you are using.

#### Windows

On Windows, you need to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use. You need a recent version of Rtools42.

#### MacOS

On macOS, you can use xcode syntax from terminal via [homebrew](https://brew.sh).

`sudo xcode-select --install`

and [Fortran v11.5](https://mac.r-project.org/tools/).

##### R package

After installing dependencies, it can install `heterogen` via remote from github .

```         
if (!require('remotes')) install.packages('remotes')
remotes::install_github("patauchi/heterogen")
```

# Introduction

`heterogen` is an R package made to create heterogeneity layers as additional information in ecological niche models.

# Getting Started

## Environmental Data

In this approach, whatever climate data can be used such as WorldClim, Chelsa, MerraClim, etc.

## Heterogeneity Layers

#### Repository

A complete set of raster layers was created to promote the use of spatial heterogeneity of environmental conditions. The datasets are available from Figshare [Link](https://doi.org/10.6084/m9.figshare.23903574.v1).

## Examples

## Warnings

#### Performance

At this time, the core of the functions are constantly modified in order to reduce the time of data processing (large dataset: 1Mx1M matrices). Stable version are available from CRAN.

## References

Atauchi, P. Joser; Peterson, A.Townsend. (202X). Incorporating spatial environmental heterogeneity as additional information on environmental and spatial context in ecological niche models. *Ecological Modelling, XXXX*.

Atauchi, P. Joser; Peterson, A. Townsend (2023). A global dataset of Heterogeneity at high resolution. figshare. Dataset. <https://doi.org/10.6084/m9.figshare.23903574.v2>
