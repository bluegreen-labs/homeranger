#' Format drivers
#'
#' Download and convert drivers for a given track.
#'
#' @param track biologging track with coordinates in lat lon (sf object)
#' @param id planetary computer id
#' @param secret planetary computer secret
#'
#' @return a nested list with the data cube or state space to use, the converted
#' X/Y coordinates relative to the data cube, and the geo-referencing information
#' for mapping and back converting the estimated locations
#' @export
#' @import rstac

hr_format_drivers <- function(track, id, secret){

  # the input is a track file with locations
  # in geographic coordinates, these define the bounding
  # box of the raster data to download

  # NOTE: add padding routine
  bbox <- sf::st_bbox(track)

  # rstac download
  # dem <- rstac_download(bbox)
  # lc <- rstac_download(bbox)
  # slope <- sf::gdal_utils()

  # resample to common resolution and combine
  # cube <- terra::rast(c(slope, slope^2, forest_cover_325m, forest_cover))

  # extract locations of gps points from the cube
  # coords <- extract(track, cube, xy = TRUE)

  # drop referencing information save externally
  #geo <- get_geo(cube)
  #cube_final <- drop_geo(cube)

  # return(cube_final, coords, geo)
}

# should contain the format drivers function
# using either the CDSE package or using the STAC data
# directly
#
# Required data is CopDEM data (10m)
# and CLMS 0.1 ha data / S2GLC data on the STAC
#
# Planetary computer provides the best support for
# data downloads it seems, as well as examples for download
#
# https://planetarycomputer.microsoft.com/docs/quickstarts/reading-stac-r/
# library(rstac)
