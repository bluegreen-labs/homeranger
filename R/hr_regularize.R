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
#' @param buffer buffer (in km, default = 20) around the bio-logging track file
#'  to accommodate for a sufficiently large home range area. Value should be
#'  adjusted to the target species.
#'
#' @returns nested list with the regularized data and observed track, also saves
#'  the used projection crs
#' @export

hr_regularize <- function(
  track,
  dem,
  lc,
  crs = "EPSG:3035",
  buffer = 20
){

  # set progression feedback to FALSE
  terra::terraOptions(progress = FALSE)

  if(missing(track)){
    stop("missing region of interest")
  }

  # calculate reprojected bounding box
  bbox <- track |>
    sf::st_transform(4326) |>
    sf::st_buffer(dist = units::as_units(buffer,"km")) |>
    sf::st_bbox()

  # download the data on the fly if data input is missing
  if(missing(dem) & missing(lc)){
    cli::cli_alert(
      "Missing matching DEM and Land Cover data - downloading data"
    )

    # download the data into the temporary directory
    hr_dl_maps(
      track,
      path = tempdir(),
      buffer = buffer,
      overwrite = TRUE
    )

    # read the data
    dem <- terra::rast(file.path(tempdir(), "DEM.tif"))
    lc <- terra::rast(file.path(tempdir(), "LC.tif"))
  }

  # re-project to equal area
  cli::cli_alert("Reproject data to CRS {crs}")
  dem <- terra::project(dem, crs, method = "bilinear")
  lc <- terra::project(lc, crs, method = "mean")

  # resample dem (30m) to lc grid (10m)
  cli::cli_alert("Resampling (DEM) data to the highest resolution (LC)")
  dem <- terra::resample(dem, lc, method = "bilinear")

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
  lc <- lc |>
    terra::crop(
      bbox
    )

  # split out mask
  m <- (lc >= 60 | lc <= 80)

  # filter out the forest data
  # and calculate tree cover density for
  # a ~310m window
  cli::cli_alert("Calculating downsampled forest cover data")
  forest <- (lc == 10)
  forest_density <- terra::focal(forest, w = 31, fun = "mean")
  forest_density_sq <- forest_density^2

  # agriculture, pasture, urban
  ag <- (lc >= 30 & lc <= 50)

  # regrowth/reforested (shrubland) (code 5322 in original data)
  shrub <- lc == 20

  # z-score conversion
  cli::cli_alert("Calculating z-score values")
  slope <- z_score(slope)
  slope_sq <- z_score(slope_sq)
  forest_density <- z_score(forest_density)
  forest_density_sq <- z_score(forest_density_sq)

  # I don't think we need z-scores for land
  # cover classes
  #ag <- z_score(ag)
  #shrub <- z_score(shrub)

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

  # return compiled dataset
  return(output)
}
