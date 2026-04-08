
# load compressed data
load(system.file("extdata/raster_maps.rda", package = "homeranger"))
load(system.file("extdata/reference_data.rda", package = "homeranger"))

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

  # run the model for these parameters
  output <- hr_predict(
    data = raster_maps,
    par = params,
    obs = system.file("extdata/Aspromonte_roedeer_traj.txt", package = "homeranger"),
    resolution = 25,
    steps = 0,
    runs = 0,
    optimization = FALSE,
    verbose = FALSE
  )

  # total residuals should add up to 0
  output$likelihood[output$likelihood == -9999] <- NA
  residuals <- sum(output$likelihood - reference_data$likelihood, na.rm = TRUE)
  print(residuals)
  expect_equal(residuals, 0)
})

test_that("test optimizations", {

  settings <- list(
    metric = hr_cost,
    control = list(
      sampler = "DEzs",
      settings = list(
        burnin = 10,
        iterations = 60
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
    drivers = raster_maps,
    obs = system.file("extdata/Aspromonte_roedeer_traj_1196.txt", package = "homeranger"),
    resolution = 25,
    settings = settings,
    parallel = FALSE
  )

  expect_type(pars, "list")
})

# test_that("test helper functions", {
# })

