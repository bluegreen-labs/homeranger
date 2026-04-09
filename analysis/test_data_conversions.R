library(terra)
library(sf)

# specify the parameters as used in the default run
# this is would be config_best_Mmem_fitting.txt
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

# read in cleaned and regularized bio-logging observations in
# lon/lat

track <-read.csv("data-raw/tracks/regularized_data_final.csv", sep = ";") |>
  na.omit() |>
  sf::st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")

# calculate bounding box around track
track_bbox <- track |>
  sf::st_transform(4326) |>
  sf::st_bbox() |>
  sf::st_as_sfc() |>
  sf::st_buffer(dist = units::as_units(20,"km"))

source("R/helpers.R")
source("R/hr_regularize.R")
source("R/hr_xy.R")
source("R/hr_dl_maps.R")

#hr_dl_maps(track, path = "analysis/")
hr_regularize(track = track, path = "analysis/test.tif")
