library(homeranger)
library(terra)

# read in the reference data, these are calculated with the
# shared original code and provide the step based likelihoods
# this output should match the output of the package for parity
reference <- read.csv("data-raw/validation/objective_function_detail.csv")
reference$likelihood[reference$likelihood == -9999] <- NA

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

# sort and subset data
# load the raster data in a data cube, and reorder the layers
# based upon the order of the coefficients in the parameter
# list - finally convert to 3D array to be passed to the
# Cpp function

original <- TRUE

if(original){
  r <- terra::rast("analysis/test_asc.tif")
  res <- terra::res(r)[1]
  r <- as.array(subset(r, names(params$coef)))
  r[is.na(r)] <- 0
  obs <- read.csv("data-raw/tracks/Aspromonte_roedeer_traj.txt") |>
    as.matrix()

} else {
  r <- terra::rast("analysis/test.tif")
  track <-read.csv("data-raw/tracks/regularized_data_final.csv", sep = ";") |>
    dplyr::rename(id = animals_id) |>
    na.omit() |>
    sf::st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")
  obs <- hr_xy(r, track)
  res <- terra::res(r)[1]
  r <- terra::as.array(r)
}

# run the model for these parameters
# in optimization mode (to check a traceable output)
# there should be ~parity as this is deterministic
output <- hr_predict(
  data = r,
  par = params,
  obs = obs,
  steps = 1,
  runs = 1,
  resolution = res,
  optimization = TRUE,
  verbose = TRUE
)

# plot the 1:1 graph - should be spot on
output$likelihood[output$likelihood == -9999] <- NA
plot(output$likelihood, reference$likelihood)
abline(0,1)

# run the model for these parameters in prediction
# mode
output <-
  hr_predict(
    data = r,
    par = params,
    obs = obs,
    resolution = res,
    steps = 10,
    runs = 2,
    verbose = TRUE
  )

# print method for hr_predict class
print(head(output$locations, 20))
plot(output)

