# a 100 runs on the latest install
#
# run time depends on the system:
# on reference system 45 - 60s

# data conversion

# r <- rast(list.files("data-raw/drivers/","*.asc", full.names = TRUE))
# data <- hr_convert_drivers(r, na_fill = 0)
#
# obs <- read.csv("data-raw/tracks/Aspromonte_roedeer_traj.txt") |>
#   #dplyr::filter(animal_id == 1196) |>
#   dplyr::mutate(across(where(is.numeric), ~dplyr::na_if(., -9999))) |>
#   dplyr::mutate(
#     x = as.integer(x / 25),
#     y = as.integer(1200 - (y / 25))
#   ) |>
#   dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(., -9999))) |>
#   as.matrix() |>
#   as.data.frame()

if(as.numeric(packageDescription("homeranger")$Version) <= 0.4){
 stop("Install the latest release or rebuild package")
}

library(homeranger)
library(terra)
library(dplyr)
library(tidyr)
source("R/hr_cost.R")

settings <- list(
  metric = "hr_cost",
  control = list(
    sampler = "DEzs",
    settings = list(
      iterations = 300
    )
  ),
  par = list(
    r_l = list(lower=10, upper=30, init = 15),
    w_l = list(lower=-1, upper=1, init = 0.5),
    r_d = list(lower=0.0001, upper=1, init = 0.5),
    w_d = list(lower=0.0001, upper=1, init = -0.5),
    r_dist = list(lower=0.0001, upper=1, init = 0.5),
    w_dist = list(lower=0.0001, upper=1, init = 0.5),
    step_length_dist = list(lower=0.0001, upper=0.1, init = 0.5),
    step_length_shape = list(lower=0.0001, upper=3, init = 1),
    threshold_approx_kernel = list(lower=6999, upper=7001, init = 7000),
    threshold_memory_kernel = list(lower=999, upper=1001, init = 1000),

    # resource selection coefficients come last
    # these are unnamed
    coef = list(
      lower = rep(-1, 6),
      upper = rep(1, 6),
      init = rep(0, 6)
    )
  )
)

load(system.file("extdata/roedeer.rda", package = "homeranger"))
load(system.file("extdata/drivers.rda", package = "homeranger"))

obs <- obs |> dplyr::filter(id == 1194)

# calibrate the model and optimize free parameters
# for only ONE individual !!
pars <- hr_fit(
    data = drivers,
    obs = obs,
    par = settings,
    parallel = FALSE,
    verbose = TRUE
)

# calibrate the model and optimize free parameters
# for only ONE individual using parallel processing across
# proposals !!
pars <- hr_fit(
  data = drivers,
  obs = obs,
  par = settings,
  parallel = TRUE,
  verbose = TRUE
)
