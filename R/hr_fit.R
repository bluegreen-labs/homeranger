#' Fit a home range model
#'
#' This is the main function that handles the
#' calibration of the home_range_cpp() C++ function
#' through hr_cost().
#'
#' @param drivers A data frame with driver data.
#' @param obs A data frame containing observational data used for model
#'  calibration.
#' @param settings A list containing model calibration settings.
#' @param optim_out A logical indicating whether the function returns the raw
#'  output of the optimization functions (defaults to TRUE).
#' @param parallel support parallel processing
#' @param ... Optional arguments passed on to the cost function specified as
#'  \code{settings$metric}.
#' @return A named list containing the calibrated parameter vector `par` and
#' the output object from the optimization `mod`. For more details on this
#' output and how to evaluate it, see \link[BayesianTools:runMCMC]{runMCMC} (also
#' \href{https://florianhartig.github.io/BayesianTools/articles/BayesianTools.html}{this post})
#' @export
#' @import BayesianTools

hr_fit <- function(
    drivers,
    obs,
    settings,
    optim_out = TRUE,
    parallel = FALSE,
    ...
){
  # predefine variables for CRAN check compliance
  lower <- upper <- out_optim <- NULL

  # check input variables
  if(missing(obs) | missing(drivers) | missing(settings)){
    stop("missing input arguments, please check all parameters")
  }

  # check data structure
  # if(is.data.frame(obs)){
  #   if (nrow(obs) == 0){
  #     warning("no validation data available, returning NA parameters")
  #     return(lapply(settings$par,
  #                   function(x) NA))
  #   }
  # }else{
  #   stop("obs must be a (nested) data.frame")
  # }

  # convert to standard cost function naming
  # cost <- eval(settings$metric)
  cost <- function(par, obs, drivers, names){
    eval(settings$metric)(
      par = par,
      obs = obs,
      data = drivers,
      names = names,
      ...
    )
  }

  # reformat parameters
  pars <- par_delist(settings$par)

  priors  <- BayesianTools::createUniformPrior(
    unlist(pars$lower),
    unlist(pars$upper),
    unlist(pars$init)
  )

  # setup the Bayes run, no message forwarding is provided
  # so wrap the function in a do.call
  setup <- BayesianTools::createBayesianSetup(
    likelihood = function(random_par) {
      do.call("cost",
              list(
                par = random_par,
                obs = obs,
                drivers = drivers,
                names = rownames(pars)
              ))
    },
    prior = priors,
    names = rownames(pars),
    parallel = parallel
  )

  # set bt control parameters
  bt_settings <- settings$control$settings

  # calculate the runs
  out <- BayesianTools::runMCMC(
    bayesianSetup = setup,
    sampler = settings$control$sampler,
    settings = bt_settings
  )

  # drop last value
  bt_par <- BayesianTools::MAP(out)$parametersMAP
  bt_par <- bt_par[1:(length(bt_par))]
  if(optim_out){
    out_optim <- list(par = bt_par, mod = out)
  }else{
    out_optim <- list(par = bt_par)
  }

  print(out_optim$par)
  #names(out_optim$par) <- rownames(pars)

  return(out_optim)
}
