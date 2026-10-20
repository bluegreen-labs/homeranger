This is the original model from the paper as listed here:
https://github.com/ranc-n/Ranc-et-al.-EcoLett-2021/tree/main

Small changes were made to make the code function, this includes adding the
chronos library to header.h and switching the Main.cpp function to its
command line setting (i.e. true).

The code is run in optimization mode primarily and returns log likelihoods.
The output of this routine is used to check for model parity in `analysis/test_predictions.R`.
The latter scripts sets the processing to optimization mode (manually) to retrieve the same
output (for a fixed set of input paramters), true predictions include a random component which
is hard to change in the original code.
