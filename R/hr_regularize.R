#' Regularize spatial input data
#'
#' Regularizes spatial input data to a common equal area grid. It also performs
#' sanity checks on the inputs to check if observations cover the spatial grid.
#' Data can be provided as arguments, or if no inputs are provided the data will
#' be downloaded in the background.
#'
#' @param dem DEM terra raster
#' @param lc land cover map (ESA world cover)
#' @param track observed locations track
#' @param crs an equal area CRS projection reference (default = EPSG:3035,
#'  or extended LAEA Europe)
#'
#' @returns nested list with the regularized data and observed track, also saves
#'  the used projection crs
#' @export

hr_regularize <- function(
    dem,
    lc,
    track,
    crs = "EPSG:3035"
){

  if(missing(track)){
    stop("missing region of interest")
  }

  # download the data on the fly if data input is missing
  if(missing(dem) & missing(lc)){
    warning("Missing matching DEM and Land Cover data - downloading data")

    # download the data into the temporary directory
    hr_dl_maps(track, path = tempdir(), overwrite = TRUE)

    # read the data
    dem <- terra::rast(file.path(tempdir(), "DEM.tif"))
    lc <- terra::rast(file.path(tempdir(), "LC.tif"))
  }

  # reproject to equal area
  dem <- terra::project(dem, crs = crs, method = "bilinear")
  lc <- terra::project(lc, crs = crs, method = "mean")

  # resample dem (30m) to lc grid (10m)
  dem <- terra::resample(dem, lc, method = "bilinear")

  # Process DEM data
  slope <- dem |>
    terra::terrain(v = "slope") |>
    terra::crop(track_bbox)
  slope_sq <- slope^2

  # Process Land Cover data
  # crop data
  lc <- lc |>
    terra::crop(track_bbox)

  # split out mask
  m <- (lc >= 60 | lc <= 80)

  # filter out the forest data
  forest <- (lc == 10)

  # agriculture
  ag <- (lc == 40 | lc == 30)

  # regrowth/reforested (shrubland)
  shrub <- lc == 20

  # output raster
  test <- c(slope, slope_sq, forest, ag, shrub)



}
