library(homeranger)

par <- list(
  metric = "hr_cost",
  control = list(
    sampler = "DEzs",
    settings = list(
      iterations = 100
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
    threshold_approx_kernel = list(lower=400, upper=402, init = 401),
    threshold_memory_kernel = list(lower=400, upper=402, init = 401),

    # resource selection coefficients come last
    # these are unnamed
    coef = list(
      lower = rep(-3, 6),
      upper = rep(3, 6),
      init = rep(0, 6)
    )
  )
)

load(system.file("extdata/roedeer.rda", package = "homeranger"))
load(system.file("extdata/drivers.rda", package = "homeranger"))

obs <-obs |>
  as.data.frame()

obs <- obs |>
  dplyr::group_split(animal_id)

# Create cluster with n cores
cl <- parallel::makeCluster(4)

# export functions / variables
parallel::clusterExport(
  cl,
  varlist = list("drivers","par", "hr_fit")
)

split_fit <- function(X){
  hr_fit(
    data = drivers,
    obs = X,
    par = par,
    parallel = FALSE
  )
}

opt_list <- parallel::parLapply(
  cl = cl,
  X = obs,
  fun = split_fit
)

parallel::stopCluster(cl)

# split out model fits
mod <- lapply(opt_list, \(x) x$mod)

## Combine the chains
out <- BayesianTools::createMcmcSamplerList(mod)

