#' Optimization cost function
#'
#' Wrapper around the low level Cpp function, avoids some of the
#' sanity checks of hr_run() thereby saving time during optimization.
#'
#' The input assumes regularized data using hr_regularize() for both
#' the spatial grid and the observed tracks.
#'
#' @param par parameters
#' @param obs reference observations
#' @param data driver matrix
#' @param names names of the parameters
#'
#' @returns log likelihood values given a data and parameter set
#' @export

hr_cost <- function(
    par,
    obs,
    data,
    names
){

  # reformat the parameters
  par_ref <- par_relist(par, names)

  # steps and runs are ignored during
  # optimization
  cost <- home_range_cpp(
    data = data,
    par = par_ref,
    trajectoryPath = obs,
    resolution = 25,
    nSimulatedSteps = 0,
    nSimulatedRuns = 0,
    optimization = TRUE,
    verbose = FALSE
  )

  # weighted mean cost can be calculated using
  # the observed values
  cost$likelihood[cost$likelihood == -9999] <- NA

  # return singe log likelihood
  return(sum(cost$likelihood, na.rm = TRUE))
}
