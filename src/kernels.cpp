// This holds the (moving window) look-up tables or kernels
// which are used to either determine the step probability
// or hold temporary local memory (these kernels are therefore
// moving windows in either space or time)

#include "header.h"

// Creates a look-up table for the movement probability
// based on the step length distribution
lookupTable iniApproxKernelStepLength (
    double distance_threshold,
    double resolution,
    double n_bins,
    double spatial_decay,
    double shape,
    double residence_p
  ){

  // create a 2D "matrix" look-up table to hold the data
  lookupTable returnValues;
  returnValues.nCells = (int) distance_threshold / resolution;
  int length = 2 * returnValues.nCells + 1;
  initialize2D(returnValues.vals, length, length);

  double dr, dc;
  double dr_sq, dc_sq;
  double dist;
  double size_bins = resolution/(2 * n_bins);
  double integral_numerator = 0;
  double integral_denominator = 0;
  double dist_correction = 0;
  double correction_factor = 0;
  double val = 0;
  double val_2 = 0;
  double sum = 0;

  for(int r = 0; r < length; r++){
		dr = (r - returnValues.nCells) * resolution;
		dr_sq = dr * dr;

		for(int c = 0; c < length; c++){
			dc = (c - returnValues.nCells) * resolution;
			dc_sq = dc * dc;
			dist = sqrt(dr_sq + dc_sq);

			if(dist > distance_threshold){
				returnValues.vals[r][c]=0;
			} else {

			// This implements the weibull distribution centered
			// on the middle of the matrix (window) for the given
			// parameters
			if(dist == 0){
  			for(int x = (-n_bins); x <= n_bins; x++){
  				for(int y = (-n_bins); y <= n_bins; y++){

  				  // distance vector for correction
  					dist_correction = pow(
  					    (pow(x * size_bins + (size_bins / 2), 2) +
  					      pow(y * size_bins + (size_bins/2), 2)),
  					  0.5);

  					if(dist_correction > 0){
  						val = exp(-1 * pow((spatial_decay * dist_correction), shape));
  						val_2 = pow((spatial_decay * dist_correction), (shape - 1));
  						integral_numerator = integral_numerator + val * val_2;
  						integral_denominator = integral_denominator + val * val_2 *
  						  (1 / dist_correction);
  						correction_factor = integral_numerator / integral_denominator;
  					}
  				}
  			}
		   dist = dist + correction_factor;
			}
			returnValues.vals[r][c] = pow((spatial_decay * dist),(shape-1)) *
			  exp(-1 * pow((spatial_decay * dist), shape)) * (1/dist);
			sum = sum + returnValues.vals[r][c];
			}
    }
  }

  // Accounting for residence probability
  for(int r = 0; r < length; r++){
  		for(int c = 0; c < length; c++){
  			returnValues.vals[r][c] = (returnValues.vals[r][c]/sum) * (1 - residence_p);
      }
  }

  returnValues.vals[returnValues.nCells + 1][returnValues.nCells + 1] =
    returnValues.vals[returnValues.nCells + 1][returnValues.nCells + 1] + residence_p;
  return returnValues;
}

// Creates a look-up table for the local calculations
// (memory is not calculated over the entire landscape but only locally)
lookupTable iniApproxKernel (
    double distance_threshold,
    double resolution,
    double spatial_decay
  ){

  // initiate 2D matrix
  lookupTable returnValues;
  returnValues.nCells=(int) distance_threshold / resolution;
  int length = 2 * returnValues.nCells + 1;
  initialize2D(returnValues.vals, length, length);

  double dr, dc;
  double dr_sq, dc_sq;
  double dist;

  for(int r = 0; r < length; r++){
		dr = (r - returnValues.nCells) * resolution;
		dr_sq = dr * dr;

		for(int c = 0; c < length; c++){
			dc = (c - returnValues.nCells) * resolution;
			dc_sq = dc * dc;
			dist = sqrt(dr_sq + dc_sq);

			if(dist > distance_threshold){
				returnValues.vals[r][c] = 0;
			}	else {
				returnValues.vals[r][c] = exp(spatial_decay * dist);
			}
    }
  }
  return returnValues;
}
