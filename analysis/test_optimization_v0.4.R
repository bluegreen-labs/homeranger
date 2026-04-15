# reference speed test for optimization of
# a 1000 runs on the reference version v0.4
# run in a terminal in a clean R session
# i.e. not in any IDE with an R project
#
# run time depends on the system:
# on reference system 45 - 60s

if(packageVersion("homeranger") != "0.4"){
  try(detach("package:homeranger", unload = TRUE))
  try(remove.packages("homeranger", lib = "~/R/x86_64-pc-linux-gnu-library/4.5"))
  remotes::install_github(
    "bluegreen-labs/homeranger@v0.4",
    upgrade = "never",
    quiet = TRUE
  )
}

library(homeranger)
library(terra)

load(system.file("extdata/raster_maps.rda", package = "homeranger"))

params <- list(
  metric = hr_cost,
  control = list(
    sampler = "DEzs",
    settings = list(
      iterations = 1000
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
# for only ONE individual!!
pars <- hr_fit(
    data = raster_maps,
    obs = "data-raw/tracks/Aspromonte_roedeer_traj_1196.txt",
    par = params,
    resolution = 25,
    parallel = FALSE
)

# plot the parameter distributions
plot(pars$mod)
