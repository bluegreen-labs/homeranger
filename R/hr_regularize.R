#' Regularize spatial input data
#'
#' Regularizes spatial input data to a common equal area grid. It also performs
#' sanity checks on the inputs to check if observations cover the spatial grid
#' and returns the
#'
#' @param data raster stack with geospatial input data
#' @param obs observed locations track
#' @param crs an equal area CRS projection reference (default = EPSG:3035,
#'  or extended LAEA Europe)
#'
#' @returns nested list with the regularized data and observed track, also saves
#'  the used projection crs
#' @export

hr_regularize <- function(
    data,
    obs,
    crs = "EPSG:3035"
){

  # reproject the data and observations to an equal area setup
  # and a common grid

  # validate input data

  # - check if all observations cover the grid layout
  # - subset layers to specified ones in parameters
  # - error when not all layers / names are there
  # - extract resolution from data geospatial information (double?)

  # subset the input data to conform to the
  # provided names of layers in the parameters
  #data_subset <- subset(data, par$names)

  # retain the extent and projection details of the underlying
  # grid setup

  # extract the cell values (XY) of the observations relative
  # to the underlying grid


}
