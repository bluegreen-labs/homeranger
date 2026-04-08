#---- validating the location extraction wrt known tracks ----
source("R/hr_xy.R")

track <-read.csv("data-raw/tracks/regularized_data_final.csv", sep = ";") |>
  rename(
    id = animals_id
  ) |>
  na.omit() |>
  sf::st_as_sf(coords = c("x", "y"), crs = "EPSG:4326")

r <- terra::rast("analysis/test.tif")

ref <- read.csv("data-raw/tracks/Aspromonte_roedeer_traj.txt") |>
  rename(
    id = animal_id
  )

ap <- read.csv(hr_xy(r, track, method = "approx"))
ex <- read.csv(hr_xy(r, track, method = "exact"))

print(head(ref |>
             dplyr::filter(
               id == 1181
             )))
print(head(ap))
print(head(ex))

s <- rast(as.array(r$slope))
plot(s)
points(ap$x/25, ap$y/25, col = "red")
points(ex$x/25, ex$y/25, col = "yellow")
points(ref$x/25, ref$y/25, col = "green")

#---- validating the data conversions to array and back -----

# internal model logic uses plain arrays / matrices the row
# order of the data is assumed to be bottom to top rather than
# top to bottom - question is how the asc files are read into
# the cpp code and if the order of the rows follows the one
# provided in the location files
m <- matrix(0, nrow = 3, ncol = 5)
m[1, ] <- 1
print(m)

r <- rast(m)
writeRaster(r, "analysis/conversion.asc", overwrite = TRUE)
r <- rast("analysis/conversion.asc")
plot(r)
points(1,1, col = "red")

print(r[1,1])
print(terra::extract(r, data.frame(1,1)))

a <- as.array(r)
print(a)
