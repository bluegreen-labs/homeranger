#' De-list parameters
#'
#' Converts nested parameter lists to a long format
#' where the vector parameter 'coef' is split in multiple
#' parts. This allows for BayesianTools priors to be generated
#' from the nested parameter files used in the model call, and
#' allows for flexible use of driver files.
#'
#' @param par home range model parameter settings
#'
#' @returns long format parameter file
#' @export

par_delist <- function(par){
  pars_no_coef <- par
  pars_no_coef$coef <- NULL
  pars_no_coef <- t(simplify2array(pars_no_coef))
  pars_coef <- simplify2array(par$coef)
  if(is.null(nrow(pars_coef))){
    pars_coef <- t(pars_coef)
    rownames(pars_coef) <- "coef_1"
  } else {
    rownames(pars_coef) <- paste("coef", 1: nrow(pars_coef),sep = "_")
  }
  as.data.frame(rbind(pars_no_coef, pars_coef))
}

#' De-list parameters
#'
#' Converts nested parameter lists to a long format
#' where the vector parameter 'coef' is split in multiple
#' parts. This allows for BayesianTools priors to be generated
#' from the nested parameter files used in the model call, and
#' allows for flexible use of driver files.
#'
#' @param par home range model parameter settings
#' @param names parameter names
#'
#' @returns long format parameter file
#' @export
par_relist <- function(par, names){
  pars_no_coef <- par[!grepl("coef", names)]
  pars_coef <- c(par[grepl("coef", names)])
  l <- c(pars_no_coef, list(pars_coef))
  names(l) <- c(names[!grepl("coef", names)],"coef")
  return(l)
}
