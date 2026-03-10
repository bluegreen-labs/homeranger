# homeranger <img src="logo.png" align="right" height="138.5"/>

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/homeranger)](https://cran.r-project.org/package=homeranger)
[![](https://cranlogs.r-pkg.org/badges/grand-total/homeranger)](https://cran.r-project.org/package=homeranger)

This package allows for the estimation and simulation of home range behaviour, including both landscape patch dynamics (step selection) and memory effects.

## How to cite this package

> Ranc and Hufkens (2026). homeranger:
> Home range modelling framework etc.
> <https://doi.org/10.5281/zenodo.XYZ>.
> and the original model description:
> Ranc, N., Cagnacci, F., Moorcroft, P.R., 2022. Memory drives the formation of animal home ranges: Evidence from a reintroduction. Ecology Letters 25, 716–728. https://doi.org/10.1111/ele.13869.

## Installation

### stable release (PENDING)

To install the current stable release use a CRAN repository:

``` r
install.packages("homeranger")
library("homeranger")
```

### development release

To install the development releases of the package run the following
commands:

``` r
if(!require(remotes)){install.packages("remotes")}
remotes::install_github("bluegreen-labs/homeranger")
library("homeranger")
```

Vignettes are not rendered by default, if you want to include additional
documentation please use:

``` r
if(!require(remotes)){install.packages("remotes")}
remotes::install_github("bluegreen-labs/homeranger", build_vignettes = TRUE)
library("homeranger")
```
## Use




## Acknowledgements

Logo background by [Mohamed Smina on Vecteezy.com](https://www.vecteezy.com/vector-art/66620012-scenic-mountain-trail-landscape-hiking-path-through-appalachian-mountains)
