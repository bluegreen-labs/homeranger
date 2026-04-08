#' Convert map / track records to col/row (X/Y) coordinates
#'
#' Converts map and track coordinates to col/row (X/Y), relative values with
#' respect to the provided map. This is an internal function to interface
#' with the C++ code which is not georeferenced internally. Note that you will
#' need to use an equal area projection for best performance. For small areas
#' lon/lat (WGS84) can be used but is not encouraged. The project should be
#' chosen based on the locale.
#'
#' @param map a \{terra\} raster map (stack)
#' @param track an \{sf\} multi-point file
#' @param method exact or approximate location conversions (default = "exact")
#'
#' @returns a nested list containing a converted map array and a list of
#'  coordinates for the provided track points relative to this array.
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
    res <- terra::res(map)

    # extract XY coordinates
    p <- terra::extract(map, track, xy = TRUE) |>
      dplyr::mutate(
        id = track$id,
        x = floor(terra::colFromX(map, .data$x) * res),
        y = floor((nrow(map) - terra::rowFromY(map, .data$y)) * res)
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

  # write to file
  filename <- file.path(tempdir(), "track.txt")
  write.csv(p, file = filename, row.names = FALSE)

  # return list of coordinates
  return(filename)
}
