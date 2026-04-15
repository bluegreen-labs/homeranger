library(terra)
library(sf)
library(homeranger)
library(dplyr)
source("R/hr_xy.R")
source("R/hr_predict.R")

r <- terra::rast("analysis/test_asc.tif")
r_row <- nrow(r)
r_col <- ncol(r)
res <- res(r)[1]

p <-read.csv("data-raw/tracks/Aspromonte_roedeer_traj.txt") |>
  dplyr::rename(id = "animal_id")

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

r <- terra::rast("analysis/test_asc.tif")
res <- terra::res(r)[1]
r <- as.array(subset(r, names(params$coef)))
r[is.na(r)] <- 0

# run the model for these parameters
# in optimization mode (to check a traceable output)
# there should be ~parity as this is deterministic
output <- hr_predict(
    data = r,
    par = params,
    obs = as.matrix(p),
    steps = 1,
    runs = 1,
    resolution = res,
    optimization = TRUE,
    verbose = TRUE
)$likelihood

# read in the reference data, these are calculated with the
# shared original code and provide the step based likelihoods
# this output should match the output of the package for parity
reference <- read.csv("data-raw/validation/objective_function_detail.csv")$likelihood
reference[reference == -9999] <- NA

# plot the 1:1 graph - should be spot on
output[output == -9999] <- NA
plot(output, reference)
abline(0,1)
