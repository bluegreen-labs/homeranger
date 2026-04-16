#' Create spatial driver data for model input
#'
#' Regularizes spatial input data to a common equal area grid. It also performs
#' sanity checks on the inputs to check if observations cover the spatial grid.
#' Data can be provided as arguments, or if no inputs are provided the data will
#' be downloaded in the background.
#'
#' It is encouraged to save the data to file, rather than use a returned object.
#' This limits the memory footprint of operations.
#'
#' @param track observed locations track
#' @param crs an equal area CRS projection reference (default = EPSG:3035,
#'  or extended LAEA Europe)
#' @param buffer buffer (in km, default = 11) around the bio-logging track file
#'  to accommodate for a sufficiently large home range area. Value should be
#'  adjusted to the target species.
#' @param window window size for resampled forest cover values, number of pixels
#'  of the window (uneven value, default = 11 or approx. 300m)
#' @param path path and filename where to save the raster map data. When using
#'  a temporary file use file.path(tempdir(), "raster_drivers.tif") or similar
#' @param na_value value with which to substitute NA values (e.g. -9999)
#' @param overwrite overwrite output if the path name is the same (default = TRUE)
#'
#' @returns nested list with the regularized data and observed track, also saves
#'  the used projection crs
#' @export

hr_drivers <- function(
  track,
  crs = "EPSG:3035",
  buffer = 7,
  window = 11,
  path,
  na_value = NA,
  overwrite = TRUE
){

  # set progression feedback to FALSE
  terra::terraOptions(progress = FALSE)

  if(missing(track)){
    stop("missing region of interest")
  }

  # calculate reprojected bounding box
  bbox <- track |>
    sf::st_transform(4326) |>
    sf::st_bbox() |>
    sf::st_as_sfc() |>
    sf::st_buffer(dist = units::as_units(buffer,"km")) |>
    sf::st_transform(crs)

  # download the data on the fly if data input is missing
  cli::cli_alert("Downloading DEM and Land Cover data")

  # download the data into the temporary directory
  hr_dl_maps(
    track,
    path = tempdir(),
    buffer = buffer + 1, # add some buffer
    overwrite = TRUE
  )

  # read the data
  dem <- terra::rast(file.path(tempdir(), "DEM.tif"))
  lc <- terra::rast(file.path(tempdir(), "LC.tif"))

  # re-project to equal area
  cli::cli_alert("Reproject data to CRS {crs}")
  dem <- terra::project(dem, crs, method = "bilinear", mask = TRUE)
  lc <- terra::project(lc, crs, method = "mean", mask = TRUE)

  # resample dem (30m) to lc grid (10m)
  # downsample
  cli::cli_alert("Resampling (DEM) data to the lowest common resolution")
  lc <- terra::resample(lc, dem, method = "modal")

  ##upsample
  #cli::cli_alert("Resampling (DEM) data to the highest resolution")
  #dem <- terra::resample(dem, lc, method = "bilinear")

  # crop to buffered track size
  cli::cli_alert("cropping data to track outline + buffer")
  dem <-  dem |> terra::crop(bbox)
  lc <-  lc |> terra::crop(bbox)

  # Process DEM data
  cli::cli_alert("Calculating slope")
  slope <- dem |>
    terra::terrain(v = "slope")
  slope_sq <- slope^2

  # Process Land Cover data
  # crop data
  cli::cli_alert("Processing land cover data")

  # filter out the forest data
  # and calculate tree cover density for
  # a ~310m window
  cli::cli_alert("Calculating downsampled forest cover data")
  forest <- (lc == 10)
  forest_density <- terra::focal(forest, w = window, fun = "mean")
  forest_density_sq <- forest_density^2

  # agriculture, pasture, urban
  ag <- (lc >= 30 & lc <= 50)

  # regrowth/reforested (shrubland) (code 5322 in original data)
  shrub <- (lc == 20)

  # z-score conversion
  cli::cli_alert("Calculating z-score values")
  slope <- z_score(slope)
  slope_sq <- z_score(slope_sq)
  forest_density <- z_score(forest_density)
  forest_density_sq <- z_score(forest_density_sq)

  # set mask
  m <- (lc >= 60 & lc <= 80)

  # mask output
  slope <- terra::mask(
    slope,
    mask = m,
    maskvalues = TRUE,
    updatevalue = NA
  )

  slope_sq <- terra::mask(
    slope_sq,
    mask = m,
    maskvalues = TRUE,
    updatevalue = NA
  )

  forest_density <- terra::mask(
    forest_density,
    mask = m,
    maskvalues = TRUE,
    updatevalue = NA
  )

  forest_density_sq <- terra::mask(
    forest_density_sq,
    mask = m,
    maskvalues = TRUE,
    updatevalue = NA
  )

  shrub <- terra::mask(
    shrub,
    mask = m,
    maskvalues = TRUE,
    updatevalue = NA
  )

  ag <- terra::mask(
    ag,
    mask = m,
    maskvalues = TRUE,
    updatevalue = NA
  )

  # output raster
  output <- c(
    slope,
    slope_sq,
    forest_density,
    forest_density_sq,
    shrub,
    ag
  ) * 1.0

  # replace NAs with another NA value
  # as required by the cpp routines
  if(!is.na(na_value)){
    output[is.na(output)] <- na_value
  }

  # retain old names
  names(output) <- c(
    "slope",
    "slope_sq",
    "tcd_325grain",
    "tcd_325grain_sq",
    "landcover_5322",
    "landcover_agri"
  )

  # return compiled dataset if no path is provided
  if(missing(path)){
    return(output)
  } else {
    terra::writeRaster(output, filename = path, overwrite = overwrite)
  }
}
