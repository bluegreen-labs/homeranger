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
source("R/hr_fit.R")

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

r <- rast(list.files("data-raw/drivers/","*.asc", full.names = TRUE))
data <- hr_convert_drivers(r, na_fill = 0)

obs <- read.csv("data-raw/tracks/Aspromonte_roedeer_traj.txt") |>
  #dplyr::filter(animal_id == 1196) |>
  dplyr::mutate(across(where(is.numeric), ~dplyr::na_if(., -9999))) |>
  dplyr::mutate(
    x = as.integer(x / 25),
    y = as.integer(1200 - (y / 25))
  ) |>
  dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(., -9999))) |>
  as.matrix()

# calibrate the model and optimize free parameters
# for only ONE individual!!
pars <- hr_fit(
    data = data,
    obs = obs,
    par = params,
    parallel = TRUE,
    verbose = TRUE
)

gc()

# calibrate the model and optimize free parameters
# for only ONE individual!!
pars <- hr_fit(
  data = data,
  obs = obs,
  par = params,
  parallel = FALSE
)
