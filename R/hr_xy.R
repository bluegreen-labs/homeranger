#' Convert map / track records to col/row (X/Y) coordinates
#'
#' Converts map and track coordinates to col/row (X/Y), relative values with
#' respect to the provided map. This is an internal function to interface
#' with the C++ code which is not georeferenced internally. Note that you will
#' need to use an equal area projection for best performance. For small areas
#' lon/lat (WGS84) can be used but is not encouraged. The project should be
#' chosen based on the locale.
#'
#' @param map a geo-referenced \{terra\} raster map (stack)
#' @param track an \{sf\} multi-point file of observations
#' @param method exact or approximate location conversions (default = "exact")
#'
#' @returns the location of the track file in the correct format (written to
#'  the temporary directory)
#' @export

hr_xy <- function(
    map,
    track,
    method = "exact"
  ){

  if(missing(map)){
    cli::cli_abort("Missing map data")
  }

  if(missing(track)){
    cli::cli_abort("Missing track data")
  }

  track <- track |>
    sf::st_transform(terra::crs(map, proj = TRUE))

  cli::cli_alert("Extracting XY coordinates from sf track object")
  if (method == "exact"){

    # extract resolution from map meta-data
    #res <- terra::res(map)

    # extract XY coordinates
    p <- terra::extract(map, track, xy = TRUE) |>
      dplyr::mutate(
        id = track$id,
        x = floor(terra::colFromX(map, .data$x)), # * res for original
        y = floor((nrow(map) - terra::rowFromY(map, .data$y)))
      ) |>
      dplyr::select(
        "id",
        "x",
        "y"
      )
  } else {
    coords <- track |>
      sf::st_coordinates() |>
      as.data.frame()

    p <- data.frame(
      id = track$id,
      x = as.integer(coords$X) - as.integer(terra::ext(map)$xmin),
      y = as.integer(coords$Y) - as.integer(terra::ext(map)$ymin)
    )
  }

  # # create index values and start and stop
  # # positions
  # p <- p |>
  #   dplyr::mutate(
  #     row_idx = 1:n()
  #   ) |>
  #   dplyr::group_by(.data$id) |>
  #   dplyr::mutate(
  #     id_start = min(row_idx) - 1, # zero indexed in C++
  #     id_end = max(row_idx) - 1
  #   )
  #
  # m <- p |> dplyr::group_by(.data$id) |>
  #   dplyr::summarize(
  #     nr_locations = n(),
  #     start = id_start[1],
  #     end = id_end[1]
  #   ) |>
  #   dplyr::arrange(
  #     .data$start
  #   )
  #
  # # create additional legacy parameters
  # # such as included map dimensions
  # output <- list(
  #   locations = data.frame(
  #     p,
  #     min_row = 0,
  #     max_row = nrow(map),
  #     min_col = 0,
  #     max_col = ncol(map)
  #   ),
  #   metrics = m
  # )

  # return output
  return(p)
}
