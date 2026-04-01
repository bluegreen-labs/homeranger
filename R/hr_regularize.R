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
  # and calculate tree cover density for
  # a ~310m window
  forest <- (lc == 10)
  forest_density <- terra::focal(forest, w = 31, fun = "mean")
  forest_density_sq <- forest_density^2

  # agriculture, pasture, urban
  ag <- (lc >= 30 & lc <= 50)

  # regrowth/reforested (shrubland) (code 5322 in original data)
  shrub <- lc == 20

  # z-score conversion
  slope <- homeranger::z_score(slope)
  slope_sq <- homeranger::z_score(slope_sq)
  forest_density <- homeranger::z_score(forest_density)
  forest_density_sq <- homeranger::z_score(forest_density_sq)
  ag <- homeranger::z_score(ag)
  shrub <- homeranger::z_score(shrub)

  # output raster
  output <- c(slope, slope_sq, forest_density, forest_density_sq, shrub, ag)

  # retain old names
  names(output) <- c(
    "slope",
    "slope_sq",
    "tcd_325grain",
    "tcd_325grain_sq",
    "landcover_5322",
    "landcover_agri"
  )

  # return terra data
  return(output)
}
