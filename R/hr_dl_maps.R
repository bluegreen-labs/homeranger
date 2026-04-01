#' Download raster maps
#'
#' This downloads the basic raster maps covering your track, including a
#' buffer zone to account for ranging behaviour. Note that this does not
#' take care of data conversions or data harmoization.
#'
#' This function does not leverage cloud optimized geotiffs via GDAL's
#' vsicurl options due to the inconsistency in the deployment across the STAC.
#' Therefore tiled data is downloaded which comes at a slight data volume
#' penalty in favour of the stability across data products.
#'
#' Please consult the Planetary computer data repository for information on the
#' STAC product collection name, the required year and the asset (subset)
#' required.
#'
#' @param track biologging track with coordinates in lat lon (sf object)
#' @param settings downloads settings, data frame specifying data
#'  products and their target year and assets on the planetary computer STAC
#' @param buffer buffer (in km, default = 20) around the bio-logging track file
#'  to accommodate for a sufficiently large home range area. Value should be
#'  adjusted to the target species.
#' @param path output path
#' @param overwrite overwrite / and re-download the data
#'
#' @return raster maps covering your area of interest as geotiff
#' @export
#' @import rstac
#' @import units

hr_dl_maps <- function(
    track,
    settings = data.frame(
      collection = c("esa-worldcover", "cop-dem-glo-30"),
      year = c(2021, 2021),
      asset = c("map", "data"),
      filename = c("LC.tif", "DEM.tif")
    ),
    buffer = 20,
    path = tempdir(),
    overwrite = TRUE
){

  # calculate bounding box around track
  track_bbox <- track |>
    sf::st_transform(4326) |>
    sf::st_bbox() |>
    sf::st_as_sfc() |>
    sf::st_buffer(dist = units::as_units(20,"km"))

  # set the stac source
  stac_source <- rstac::stac(
    "https://planetarycomputer.microsoft.com/api/stac/v1"
  )

  # list the whole collection
  collections_query <-
    stac_source |>
    rstac::collections()

  # list all available data collections
  available_collections <- lapply(
    rstac::get_request(collections_query)$collections,
    \(x) x$id
  ) |>
    do.call(what = rbind)

  # check if both requested collations are in the STAC
  if(!all(settings$collections %in% available_collections)){
    stop("wrong collection name")
  }

  # cycling through all line entries of the products
  # in the settings data frame to download required products
  apply(settings, 1, function(s){

    # temporary directory
    tmp_dir <- file.path(tempdir(), "rstac", s["collection"])
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

    # check if file exists, if so skip
    filename <- file.path(path, s["filename"])

    if(file.exists(filename) & !overwrite){
      warning(sprintf("file %s exists, skipping", filename))
      return(NULL)
    }

    # query all available data from 2000 until now
    # this should cover most required products which are most
    # often single static layers
    stac_query <- rstac::stac_search(
      q = stac_source,
      collections = s["collection"],
      datetime = sprintf("%s-01-01/%s-12-31", s["year"], s["year"]),
      bbox = track_bbox |> sf::st_bbox(),
      limit = 999
    ) |> rstac::get_request()

    signed_stac_query <- rstac::items_sign(
      stac_query,
      rstac::sign_planetary_computer()
    )

    rstac::assets_download(
      signed_stac_query,
      s["asset"],
      output_dir = tmp_dir,
      overwrite = TRUE,
      progress = TRUE
    )

    # list downloaded files
    image_files <- list.files(
      tmp_dir,
      "*.tif",
      recursive = TRUE,
      full.names = TRUE
    )

    # merge data if the returned data
    # contains more than one file
    if(length(image_files) > 1){
      r_l <- lapply(image_files, terra::rast) |> terra::sprc()
      r <- terra::merge(r_l)
    } else {
      # read in the image file
      r <- terra::rast(image_files)
    }

    # crop data to size
    r <- terra::crop(r, track_bbox)

    # write to file, ignore warnings
    suppressWarnings(
      terra::writeRaster(
        r,
        filename = filename,
        overwrite = overwrite
      )
    )

    # delete temporary files
    #file.remove(image_files)
    }
  )

  # invisible output
  invisible()
}
