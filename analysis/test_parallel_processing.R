library(homeranger)

settings <- list(
  metric = "hr_cost",
  control = list(
    sampler = "DEzs",
    settings = list(
      iterations = 300
    )
  ),
  par = list(
    r_l = list(lower=10, upper=30, init = 15),
    w_l = list(lower=-1, upper=1, init = 0.5),
    r_d = list(lower=0.0001, upper=1, init = 0.5),
    w_d = list(lower=0.0001, upper=1, init = -0.5),
    r_dist = list(lower=0.0001, upper=1, init = 0.5),
    w_dist = list(lower=0.0001, upper=1, init = 0.5),
    step_length_dist = list(lower=0.0001, upper=0.1, init = 0.5),
    step_length_shape = list(lower=0.0001, upper=3, init = 1),
    threshold_approx_kernel = list(lower=6999, upper=7001, init = 7000),
    threshold_memory_kernel = list(lower=999, upper=1001, init = 1000),

    # resource selection coefficients come last
    # these are unnamed
    coef = list(
      lower = rep(-1, 6),
      upper = rep(1, 6),
      init = rep(0, 6)
    )
  )
)

load(system.file("extdata/roedeer.rda", package = "homeranger"))
load(system.file("extdata/drivers.rda", package = "homeranger"))

obs <-obs |>
  as.data.frame()

obs <- obs |>
  dplyr::group_split(id)

# Create cluster with n cores
cl <- parallel::makeCluster(4)

# export functions / variables
parallel::clusterExport(
  cl,
  varlist = list("drivers","settings", "hr_fit")
)

split_fit <- function(X){
  hr_fit(
    data = drivers,
    obs = X,
    par = settings,
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

