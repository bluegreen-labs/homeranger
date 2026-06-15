# a 100 runs on the latest install
#
# run time depends on the system:
# on reference system 45 - 60s

if(as.numeric(packageDescription("homeranger")$Version) <= 0.4){
 stop("Install the latest release or rebuild package")
}

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

# optimize the parameters
pars <- hr_fit(
    data = drivers,
    obs = obs,
    par = settings,
    parallel = TRUE,
    verbose = TRUE
)

# save compressed parameters
saveRDS(pars, "analysis/parameters_full_run.rds", compress = "xz")
