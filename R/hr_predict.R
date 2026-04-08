#' Title
#'
#' @param par parameters
#' @param data data
#' @param obs observed track
#' @param resolution resolution
#' @param steps steps
#' @param runs runs
#' @param optimization optimize
#' @param verbose feedback
#'
#' @returns simulated track
#' @export

hr_predict <- function(
  data,
  par,
  obs,
  resolution,
  steps = 1,
  runs = 1,
  optimization = FALSE,
  verbose = TRUE
){

  # remove coef names (only needed for regularization)
  # saves time as complex Rcpp/arma vector types are time expensive
  par$names <- NULL

  # call low level cpp model function
  output <- home_range_cpp(
    data = data,
    par = par,
    trajectoryPath = obs,
    resolution = resolution,
    nSimulatedSteps = steps,
    nSimulatedRuns = runs,
    optimization = optimization,
    verbose = verbose
  )

  # if returning the log likelihood from an optimization
  # run return this early
  if (!optimization){

    # set -9999 values to NA
    output$locations[output$locations == -9999] <- NA
    output$locations$id <- as.factor(output$locations$id)

    # class assignment
    class(output) <- c("hr_predict")
    return(output)

  } else {
    # set -9999 values to NA
    output$likelihood[output$likelihood == -9999] <- NA
    return(output)
  }
}
