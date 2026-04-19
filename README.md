# homeranger <img src="logo.png" alt = "logo" align="right" height="138" width = "138">

This package allows for the estimation and simulation of home range behaviour, including both landscape patch dynamics (step selection) and memory effects. The mechanistic movement model allows you to: "(1) quantify the role of memory in the movements of a large mammal reintroduced into a novel environment, and (2) predict observed patterns of home range emergence in this experimental setting" Ranc et al. (2022).

## How to cite this package

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

The underlying idea of the homeranger package and the underlying mechanistic
movement model is to approximate the behaviour of observed individuals using the 
provided model, an underlying spatial grid of driver data (i.e. a landscape
through which the individual moves) and a set of optimized parameters (knobs 
to tweak) to influence the underlying model.

### model fitting

The first step is therefore often optimizing the parameters of the mechanistic
movement model by calculating a best fit to observed values Running a model fit
requires some initial data pre-processing and has to adhere to the
required formats, mainly:

 - georeferenced (equal area) raster driver data
 - georeferenced observation tracks for individuals (can be multiples)

#### data pre-processing

Ideally, driver data is geo-referenced data downloaded through the download
routines or otherwise constructed (see advanced vignette). This data combined
with the parameter file can be passed to the `hr_convert_drivers()` function
to convert the data to the desired input format.

```r
# load georeferenced data as a {terra} SpatRaster file
r <- terra::rast("your_input_raster")

# convert the raster data combined with the parameters
# to the desired format
data <- hr_convert_drivers(r)
```

Geo-referenced observation data (an {sf} object) can be correctly oriented within 
the raw data array using the `hr_xy()` function, which extracts the rows and
columns within the data array.

```r
# reading in a text file with coordinate columns x_col and y_col
# in lon/lat format, make sure that the only other remaining 
# column is the individuals ID
track <- read.csv('file.csv') |>
  sf::st_as_sf(
    coords=c("x_col", "y_col"),
    crs=4326
  )

# combined with the driver raster data above extract the rows and columns
# relative to this raster data
obs <- hr_xy(r, track)
```

The final step in the parameter optimization is fitting the observed data to
the model output. This is an iterative process. Therefore we need to specify
upper and lower bounds for the parameters (which have to reflect reasonable
physical properties if these can be easily interpreted).

```r
# set optimization parameter ranges as well as the control parameters
# for the optimization algorithm
params <- list(
  metric = "hr_cost",
  control = list(
    sampler = "DEzs",
    settings = list(
      iterations = 7 * 10
    )
  ),
  par = list(
    r_l = list(lower=0.0001, upper=1, init = 0.5),
    w_l = list(lower=0.0001, upper=1, init = 0.5),
    r_d = list(lower=0.0001, upper=1, init = 0.5),
    w_d = list(lower=-1, upper=-0.0001, init = -0.5),
    r_dist = list(lower=0.0001, upper=1, init = 0.5),
    w_dist = list(lower=0.0001, upper=1, init = 0.5),
    step_length_dist = list(lower=0.0001, upper=0.1, init = 0.5),
    step_length_shape = list(lower=0.3, upper=3, init = 1),
    threshold_approx_kernel = list(lower=300, upper=302, init = 301),
    threshold_memory_kernel = list(lower=300, upper=302, init = 301),

    # resource selection coefficients come last
    # these are unnamed
    coef = list(
      lower = rep(-3, 6),
      upper = rep(3, 6),
      init = rep(0, 6)
    )
  )
)

# calibrate the model parameters
pars <- hr_fit(
    data = data,
    obs = obs,
    par = params,
    parallel = FALSE,
    verbose = TRUE
)

# plot the model parameter distributions
plot(pars$mod)
```

Note that the option for parallel processing exists, but only provides
significant speed gains when optimizing model parameters for multiple
individuals at the same time or very long observation tracks.

### model prediction

The data pre-processing follows the same steps as the model fitting above. The
model parameters are either provided by the fitting routine above (in this case
stored in `pars$par`).

```r
# the model parameters in the pars$par have the following structure
params <- list(
  r_l = 27.5332236990522,
  w_l = 0,
  r_d = 0.018257876686841,
  w_d = 0.9999,
  r_dist = 0.0412536435305482,
  w_dist = 0.9999,
  step_length_dist = 0.00216275705935606,
  step_length_shape = 1.14267311221975,
  threshold_approx_kernel = 7000,
  threshold_memory_kernel = 1000,

  # resource selection coefficients should be
  # a named list for driver data layer validation
  # and correct data processing
  coef = c(
    "slope" = 0.272835968106296,
    "slope_sq" = -0.093687792157105,
    "tcd_325grain"= 0.177991482087775,
    "tcd_325grain_sq" = -0.140639949444926,
    "landcover_5322" = 0.591063382485486,
    "landcover_agri" = -0.811974081226742
  )
)

# simulate 100 steps from the starting position and this for two randomizations
output <- hr_predict(
  data = data,
  par = params,
  obs = obs,
  steps = 100,
  runs = 2
)

# plot the model simulation output using
plot(output)
```

## Acknowledgements

The original model development was supported by the Harvard University Graduate Fellowship and a Fondazione Edmund Mach International Doctoral Programme Fellowship, with data and code made available at the Zenodo Digital Repository (https://doi.org/10.5281/zenodo.5189835 and https://doi.org/10.5281/zenodo.5208215, repository). The logo background was created by [Mohamed Smina and freely shared through Vecteezy.com](https://www.vecteezy.com/vector-art/66620012-scenic-mountain-trail-landscape-hiking-path-through-appalachian-mountains). Refactoring and further model development was supported by the French Office for Biodiversity on the TRANSLOC project (OF-25-0032).
