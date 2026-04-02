library(terra)
library(sf)

# read in cleaned and regularized bio-logging observations in
# lon/lat

track <-read.csv("data-raw/tracks/regularized_data_final.csv", sep = ";") |>
  na.omit() |>
  sf::st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")

# calculate bounding box around track
track_bbox <- track |>
  sf::st_transform(4326) |>
  sf::st_bbox() |>
  sf::st_as_sfc() |>
  sf::st_buffer(dist = units::as_units(20,"km"))

source("R/helpers.R")
source("R/hr_regularize.R")
source("R/hr_dl_maps.R")

test <- hr_regularize(track = track_bbox)
