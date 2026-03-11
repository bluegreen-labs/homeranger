
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

  # load compressed data
  load(system.file("extdata/raster_maps.rda", package = "homeranger"))
  load(system.file("extdata/reference_data.rda", package = "homeranger"))

  # run the model for these parameters
  output <- hr_predict(
    data = raster_maps,
    par = params,
    obs = system.file("extdata/Aspromonte_roedeer_traj.txt", package = "homeranger"),
    resolution = 25,
    optimization = FALSE,
    verbose = FALSE
  )

  # total residuals should add up to 0
  output$likelihood[output$likelihood == -9999] <- NA
  residuals <- sum(output$likelihood - reference_data$likelihood, na.rm = TRUE)
  print(residuals)
  expect_equal(residuals, 0)
})

# test_that("test optimizations", {
# })
#
# test_that("test helper functions", {
# })

