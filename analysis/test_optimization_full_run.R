# Recreating the optimization of the Ranc et al. paper.

# Testing the stability of long runs and large optimization
# jobs, this run is internally parallelized using Bayesiantools
# and uses a maximum of three (3) threads.

library(homeranger)
library(terra)
library(dplyr)
library(tidyr)

load(system.file("extdata/roedeer.rda", package = "homeranger"))
load(system.file("extdata/drivers.rda", package = "homeranger"))

settings <- list(
  metric = "hr_cost",
  control = list(
    sampler = "DEzs",
    settings = list(
      iterations = 80000
    )
  ),
  par = list(
    r_l = list(lower=0.0001, upper=1, init = 0),
    w_l = list(lower=0.0001, upper=1, init = 0),
    r_d = list(lower=0.0001, upper=1, init = 0.1),
    w_d = list(lower=-1, upper=1, init = -0.0001),
    r_dist = list(lower=0.0001, upper=1, init = 0.1),
    w_dist = list(lower=0.0001, upper=1, init = 0.1),
    step_length_dist = list(lower=0.0001, upper=0.1, init = 0.5),
    step_length_shape = list(lower=0.3, upper=3, init = 1),
    threshold_approx_kernel = list(lower=6999, upper=7001, init = 7000),
    threshold_memory_kernel = list(lower=999, upper=1001, init = 1000),

    # resource selection coefficients come last
    # these are unnamed
    coef = list(
      lower = rep(-3, 6),
      upper = rep(3, 6),
      init = rep(0, 6)
    )
  )
)

# optimize the parameters
pars <- hr_fit(
    data = drivers,
    obs = obs,
    par = settings,
    parallel = TRUE,
    verbose = TRUE
)

plot(pars$mod)

# save compressed parameters
#saveRDS(pars, "analysis/parameters_full_run.rds", compress = "xz")
