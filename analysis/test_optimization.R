library(homeranger)
load(system.file("extdata/raster_maps.rda", package = "homeranger"))
source("R/hr_fit.R")

params <- list(
  metric = hr_cost,
  control = list(
    sampler = "DEzs",
    settings = list(
      burnin = 3,
      iterations = 4
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
    threshold_approx_kernel = list(lower=3000, upper=10000, init = 7000),
    threshold_memory_kernel = list(lower=3000, upper=10000, init = 1000),

    # resource selection coefficients come last
    # these are unnamed
    coef = list(
      lower = rep(-3, 6),
      upper = rep(3, 6),
      init = rep(0, 6)
    )
  )
)

# calibrate the model and optimize free parameters
pars <- hr_fit(
  data = raster_maps,
  obs = system.file("extdata/Aspromonte_roedeer_traj_1196.txt", package = "homeranger"),
  par = params,
  resolution = 25,
  parallel = FALSE
)

# plot the parameter distributions
plot(pars$mod)
