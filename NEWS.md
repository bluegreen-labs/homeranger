## homeranger v0.5

- traps out of range memory issue around L261, which breaks optimization
- implemented proper (internal) parallelization
- included new functions for data formatting `hr_convert_drivers()`
- new data formatting allows for the removal of the `resolution` parameter in the `hr_fit()` and `hr_predict()` functions.
- updated vignettes to include data formatting and parallelization instructions (internal + external)
- fixed final potential memory leak issue

## homeranger v0.4

- raster driver regularization (downloads)
- xy extraction based on an {sf} point object
- consolidation function parameters + docs

## homeranger v0.3

- integrating both kernel estimation and simulation modes
- setting an R compatible random number generator
- adding additional benchmark functionality (not integrated into package formally only on github)

## homeranger v0.2

- adding data download function `hr_raster_maps()`
- adding documentation for data formatting
- adding a coordinate conversion function `hr_xy()`

## homeranger v0.1 (alpha)

This is the first iteration of the {homeranger} package. A framework for home range estimation and modelling. This release contains the refactored version of the original C++ code into an R package (with a C++ back-end, through Rcpp).
