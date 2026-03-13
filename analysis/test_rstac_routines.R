# testing rstac Planetary computer downloads

library("sf")
library("rstac")
library("terra")
source("R/hr_raster_maps.R")

ashe <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))[1, ]
hr_raster_maps(track = ashe, path = "analysis/")
