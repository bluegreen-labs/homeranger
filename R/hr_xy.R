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
#' @param crs a valid CRS description (default = "+proj=longlat +datum=WGS84")
#'
#' @returns a nested list containing a converted map array and a list of
#'  coordinates for the provided track points relative to this array.
#' @export

hr_xy <- function(
    map,
    track = NULL,
    crs = "+proj=longlat +datum=WGS84"
  ){

  # reproject to equal area locations
  map <- terra::project(map, crs)

  if(!is.null(track)){
    track <- track |>
      sf::st_transform(crs) |>
      terra::vect()

    # extract XY coordinates
    p <- terra::extract(map, track, xy = TRUE) |>
      dplyr::mutate(
        col = terra::colFromX(map, .data$x),
        row = terra::rowFromY(map, .data$y)
      ) |>
      dplyr::select(
        "col", "row"
      )
  } else {
    p <- NULL
  }

  # convert to array
  r <- as.array(r)
  r[is.na(r)] <- 0

  # return nested list of data
  list(raster = r, coordinates = p)
}
