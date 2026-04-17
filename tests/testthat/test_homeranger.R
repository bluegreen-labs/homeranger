
# load compressed data
load(system.file("extdata/raster_maps.rda", package = "homeranger"))
load(system.file("extdata/reference_data.rda", package = "homeranger"))

library(tidyr)
library(terra)
library(homeranger)

test_that("validate model run", {
  # set parameters
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
    coef = c(
      "slope" = 0.272835968106296,
      "slope_sq" = -0.093687792157105,
      "tcd_325grain"= 0.177991482087775,
      "tcd_325grain_sq" = -0.140639949444926,
      "landcover_5322" = 0.591063382485486,
      "landcover_agri" = -0.811974081226742
    )
  )

  # create desired input format
  data <- list(data = raster_maps, resolution = 25)

  # read in data and convert to matrix
  obs <- read.csv(
    system.file("extdata/Aspromonte_roedeer_traj.txt", package = "homeranger")) |>
    dplyr::mutate(across(where(is.numeric), ~dplyr::na_if(., -9999))) |>
    dplyr::mutate(
      x = as.integer(x / 25),
      y = as.integer(1200 - (y / 25))
    ) |>
    dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(., -9999))) |>
    as.matrix()

  # run the model for these parameters
  output <- hr_predict(
    data = data,
    par = params,
    obs = obs,
    optimization = TRUE, # returns the log-likelihood values for the given data
    verbose = FALSE
  )

  # total residuals should add up to 0
  output$likelihood[output$likelihood == -9999] <- NA
  residuals <- sum(output$likelihood - reference_data$likelihood, na.rm = TRUE)

  # include tolerance
  residuals <- ifelse(residuals < 0.00001, 0, 1)
  expect_equal(residuals, 0)
})

test_that("test optimizations", {

  params <- list(
    metric = "hr_cost",
    control = list(
      sampler = "DEzs",
      settings = list(
        burnin = 3,
        iterations = 10
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

  # create desired input format
  data <- list(data = raster_maps, resolution = 25)

  obs <- read.csv(
    system.file("extdata/Aspromonte_roedeer_traj.txt", package = "homeranger")) |>
    dplyr::filter(animal_id == 1196) |>
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
    parallel = FALSE
  )
  expect_type(pars, "list")

  pars_par <- hr_fit(
    data = data,
    obs = obs,
    par = params,
    parallel = TRUE
  )
  expect_type(pars_par, "list")
})

# test_that("test helper functions", {
# })

