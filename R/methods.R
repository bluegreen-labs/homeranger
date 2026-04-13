#' Plotting method for prediction results
#'
#' @param x input data generated using the pr_fit() function
#' @param ... additional parameters to pass
#' @return a plot (map) of predicted statistics
#' @keywords model, behaviour
#' @export
#'
#' @importFrom rlang .data

# plot data if requested
plot.hr_predict <- function(x, ...){
  stopifnot(inherits(x, "hr_predict"))

  map_data <- terra::rast(x$resources[[1]])
  locations <-x$locations
  locations_sf <- locations |>
    stats::na.omit() |>
    sf::st_as_sf(coords = c("col", "row"))

  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(
      data = map_data
    ) +
    ggplot2::geom_sf(
      data = locations_sf,
      ggplot2::aes(
        colour = .data$id
      )
    ) +
    ggplot2::facet_wrap(~ .data$run)

  # output plot
  print(p)
}
