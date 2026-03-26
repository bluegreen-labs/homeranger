# This code runs the original model, using the original files and
# settings. Only data paths are altered. No seed value was provided
# so simulations will include a degree of randomness, but general
# patterns should correspond.
#
# This setup is only tested on Linux.

# create a benchmarking folder to compile the code in, not to mess with
# the file structure or upload binaries.
bin_path <- file.path(tempdir(),"roedeer")
dir.create(bin_path, recursive = TRUE)

# copy all code to the bin_path folder
files <- list.files("benchmark/src/", "*", full.names = TRUE)
file.copy(files, bin_path, overwrite = TRUE)

# set config to use
config_file <- "./data-raw/config/config_best_Mmem_fitting.txt"
config <- read.delim(config_file, sep = ";", header = TRUE)

# set parameters (set trailing /)
config$value[1] <- here::here("benchmark/")
config$value[2] <- here::here("data-raw/drivers/")
config$value[3] <- here::here("data-raw/tracks/Aspromonte_roedeer_traj.txt")

# save new config to tmp location
write.table(
  config,
  file = file.path(bin_path, "config.txt"),
  sep = ";",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# compile the binary in bin_path
system(
  sprintf(
   "cd %s; g++ -std=c++11 -O3 Configuration.cpp Dist_lookup_table.cpp \\
   Import_traj.cpp Launch_arena.cpp Likelihood.cpp Main.cpp Patch_dynamics.cpp \\
   Writing_outputs.cpp -o redistribution_kernel_roedeer",
   bin_path
  )
)

# execute the model using the adjusted configuration files
system(
  sprintf(
    "cd %s; ./redistribution_kernel_roedeer -config %s/config.txt",
    bin_path,
    bin_path
  )
)
