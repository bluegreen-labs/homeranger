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
