#' Converts data to a consistent format
#'
#' Data conversion tool to generate a consistent
#' input format specific to the package
#'
#' @param r a terra SpatRaster
#' @param par a parameter file
#' @param na_fill fill value for NA values (default = 0)
#'
#' @returns a nested list of data to be used in model prediction and fitting
#' @export

hr_convert_drivers <- function(
    r,
    par,
    na_fill = 0
){

  if(!missing(par)){
    cli::cli_alert_info("Sorting and subsetting raster layer according to parameters")
    r <- terra::subset(r, names(par$coef))
  }

  # extract resolution
  resolution <- terra::res(r)[1]

  # list of layer names (in original order)
  layers <- terra::names(r)

  # extract projection meta-data
  crs <- terra::crs(r)

  # convert to array
  # and convert the NA values
  # to the desired fill value
  r <- terra::as.array(r)
  r[is.na(r)] <- na_fill

  # return everything as a nested list
  return(
    list(
      resolution = resolution,
      crs = crs,
      layer = layers,
      data = r
    )
  )
}
