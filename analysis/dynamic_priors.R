settings <- list(
  metric = hr_cost,
  control = list(
    maxit = 1
  ),
  par = list(
    r_l = list(lower=0, upper=100, init = 27),
    w_l = list(lower=0, upper=100, init = 27),
    r_d = list(lower=0, upper=1, init = 0.5),
    w_d = list(lower=0, upper=1, init = 0.5),
    r_dist = list(lower=0, upper=1, init = 0.5),
    w_dist = list(lower=0, upper=1, init = 0.5),
    step_length_dist = list(lower=0, upper=1, init = 0.5),
    step_length_shape = list(lower=0, upper=3, init = 1),
    threshold_approx_kernel = list(lower=0, upper=10000, init = 5000),
    threshold_memory_kernel = list(lower=0, upper=10000, init = 1000),
    coef = list(
      lower = rep(-1, 6),
      upper = rep(1, 6),
      init = rep(0, 6)
    )
  )
)

# reformat parameters
pars <- as.data.frame(do.call("rbind", settings$par), row.names = FALSE)

priors  <- BayesianTools::createUniformPrior(
  unlist(pars$lower),
  unlist(pars$upper),
  unlist(pars$init)
)
