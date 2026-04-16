//' @name home_range_cpp
//' @title Main home range estimation function
//' @description C++ back-end function doing all the heavy lifting
//' @param data data cube of the state space
//' @param par parameters for model execution or fitting
//' @param trajectoryPath location of the track file
//' @param resolution resolution (dimensions) of the state space (as not georeferenced)
//' @param nSimulatedSteps number of simulated steps to run
//' @param nSimulatedRuns number of runs to simulate
//' @param optimization return optimization outputs (TRUE or FALSE)
//' @param verbose provide verbose output (TRUE or FALSE)
//' @export
//'
//' @return simulated model locations or a target value for optimization
//'  (cost function)

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "header.h"
#include <chrono>

using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
List home_range_cpp(
    arma::cube data,
    List par,
    arma::mat locations,
    int resolution,
    int nSimulatedSteps,
    int nSimulatedRuns,
    bool optimization,
    bool verbose
  )
{

  // Environment setup ---------------------------------------------------------
  time_t iniTime = time(0);

  // define local variables
  double familiarity, sum_weights;
  int ind, ite, rowCount;
  int currentCol,currentRow,nextCol, nextRow, releaseCol, releaseRow;
  int minR, maxR, minC, maxC;
  int minRmem, maxRmem, minCmem, maxCmem;
  int lagR, lagC, lagRmem, lagCmem;
  int looktableR, looktableC;
  double weightR, weightW;
  int nRows, nCols;

  // get arena rows and columns
  // (pixels), together with the
  // resolution (m/pixel) this gives
  // the absolute scale of the scene
  nRows = data.n_rows;
  nCols = data.n_cols;

  // Read in the parameters and arguments
  // from the R list passed as an argument

  // step selection parameters
  // the use of arma::vec or Rcpp::NumericVector allows
  // for flexible parameter assignments without additional
  // logic checking overhead. Note that these are vectors
  // and take positional arguments. If only one value is
  // expected you need to query the [0] position which by
  // default is a double unless otherwise specified
  arma::vec selectionCoef = par["coef"];
  arma::vec stepLengthDist = par["step_length_dist"];
  arma::vec stepLengthShape = par["step_length_shape"];

  // kernel settings
  arma::vec thresholdApproxKernel = par["threshold_approx_kernel"];
  arma::vec thresholdMemoryKernel = par["threshold_memory_kernel"];
  arma::vec memoryRL = par["r_l"];
  arma::vec memoryWL = par["w_l"];

  arma::vec memoryRD = par["r_d"];
  double memoryRD_cplm = (1.0 - memoryRD[0]);
  arma::vec memoryWD = par["w_d"];
  double memoryWD_cplm = (1.0 - memoryWD[0]);

  arma::vec memoryRDist = par["r_dist"];
  arma::vec memoryWDist = par["w_dist"];

  // might want to just flip the sign
  // on the parameters instead!!
  memoryRDist = memoryRDist * (-1.0);
  memoryWDist = memoryWDist * (-1.0);

  // initiate arma matrices based on the size
  // of the input drivers
  auto begin = std::chrono::high_resolution_clock::now();

  // initiate original structure with nested vector layers
  ArraysDynamics Arena = launchArena(
    nRows,
    nCols
  );

  // Initialize the result matrix
  arma::mat arrayResourceSelection = exp(data.slice(0) * selectionCoef(0));

  // Iterate through the remaining slices (starting from the second slice, index 1).
  // For each slice, scale it by the corresponding vector element and
  // then perform an element-wise multiplication
  for (arma::uword i = 1; i < data.n_slices; ++i) {
    arma::mat slice = exp(data.slice(i) * selectionCoef(i));
    arrayResourceSelection %= slice;
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto dur = end - begin;
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();

  if(verbose){
    Rcout << "Reading data table, run time: "<< ms << " ms" << std::endl;
  }

  // deep copy
  NumericMatrix res = wrap(arrayResourceSelection);

  // convert from arma matrix
  for (int r = 0; r < nRows; ++r) {
    for (int c = 0; c < nCols; ++c) {
      Arena.arrayResourceSelection[r][c] = res(r, c);
    };
  };

  begin = std::chrono::high_resolution_clock::now();

  lookupTable stepLengthKernel = iniApproxKernelStepLength(
    thresholdApproxKernel[0],
    resolution,
    500,
    stepLengthDist[0],
    stepLengthShape[0],
    0
  );

  lookupTable r_memoryKernel = iniApproxKernel(
    thresholdMemoryKernel[0],
    resolution,
    memoryRDist[0]
  );

  lookupTable w_memoryKernel = iniApproxKernel(
    thresholdMemoryKernel[0],
    resolution,
    memoryWDist[0]
  );

  end = std::chrono::high_resolution_clock::now();
  dur = end - begin;
  ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
  if(verbose){
    Rcout << "lookup tables, run time: "<< ms << " ms" << std::endl;
  }

  // loading trajectory files
  structTrajectory Traj = launchTrajectoryCoordinatesMat(
    locations,
    0,
    0,
    nRows,
    nCols,
    r_memoryKernel.nCells
  );

  structSummaryTraj TrajMetrics = getTrajectoryMetrics(Traj);
  int nAnimals = TrajMetrics.animalId.size();

  // Start model run -----------------------------------------------------------
  if(verbose){
    begin = std::chrono::high_resolution_clock::now();
    if(optimization){
      Rcout << "Optimization mode" << std::endl;
    } else {
      Rcout << "Simulation mode" << std::endl;
    }
  }

  // initiate memory for the simulated trajectory
  // make conditional?
  structTrajectorySimul totalTraj;
  totalTraj.animalId.reserve(TrajMetrics.totalLength+nAnimals*nSimulatedSteps*nSimulatedRuns);
  totalTraj.run.reserve(TrajMetrics.totalLength+nAnimals*nSimulatedSteps*nSimulatedRuns);
  totalTraj.col.reserve(TrajMetrics.totalLength+nAnimals*nSimulatedSteps*nSimulatedRuns);
  totalTraj.row.reserve(TrajMetrics.totalLength+nAnimals*nSimulatedSteps*nSimulatedRuns);

  double sumAttractionVectors,randNumber;
  int totalCount;
  int rstart, rend, cstart, cend;
  totalCount = -1;
  rowCount=-1;
  lagR=0;
  lagC=0;
  lagRmem=0;
  lagCmem=0;

  // START OF INDIVIDUAL LOOP
  for(ind=0;ind<nAnimals;ind++){

    // Re-initialization of the arena
    // clearing reference and working memory
    arena_renewal(Arena, 0, 0);

    // START OF RELOCATION ITERATIONS
    // RELOCATION LOOP (from release to 2nd from last point)
    for(ite = 0; ite < TrajMetrics.individualCount[ind] - 1; ite++){

      rowCount = rowCount + 1;
      currentCol = Traj.col[rowCount];
      currentRow = Traj.row[rowCount];

      nextCol = Traj.col[rowCount + 1];
      nextRow = Traj.row[rowCount + 1];

      if(optimization){
        minR = currentRow - stepLengthKernel.nCells;
        maxR = currentRow + stepLengthKernel.nCells + 1;
        minC = currentCol - stepLengthKernel.nCells;
        maxC = currentCol + stepLengthKernel.nCells + 1;

        // check boundary conditions
        if(minR < 0){lagR = minR; minR = 0;}
        if(maxR > Arena.nRows){maxR = Arena.nRows;}
        if(minC < 0){lagC = minC; minC = 0;}
        if(maxC > Arena.nCols){maxC = Arena.nCols;}
      }

      minRmem = currentRow - r_memoryKernel.nCells;
      maxRmem = currentRow + r_memoryKernel.nCells + 1;
      minCmem = currentCol - r_memoryKernel.nCells;
      maxCmem = currentCol + r_memoryKernel.nCells + 1;

      // check boundary conditions
      lagRmem = 0;
      lagCmem = 0;
      if(minRmem < 0){lagRmem = minRmem; minRmem = 0;}
      if(maxRmem > Arena.nRows){maxRmem = Arena.nRows;}
      if(minCmem < 0){lagCmem = minCmem; minCmem = 0;}
      if(maxCmem > Arena.nCols){maxCmem = Arena.nCols;}

      // 1. Memory dynamics
      //"full" memory decay
      rstart = Traj.minRowMem[rowCount];
      rend = Traj.maxRowMem[rowCount];
      cstart = Traj.minColMem[rowCount];
      cend = Traj.maxColMem[rowCount];

      if(rstart < 0){rstart = 0;}
      if(rend > Arena.nRows){rend = Arena.nRows;}
      if(cstart < 0){cstart = 0;}
      if(cend > Arena.nCols){cend = Arena.nCols;}

      for(int r = rstart; r < rend; r++){
        for(int c = cstart; c < cend; c++){
          Arena.arrayMemoriesRef[r][c] = Arena.arrayMemoriesRef[r][c] * memoryRD_cplm;
          Arena.arrayMemoriesWork[r][c] = Arena.arrayMemoriesWork[r][c] * memoryWD_cplm;
        }
      }

      // If current timestep has valid coordinates
      if(currentCol != -9999){

        // At cells within neighborhood...
        for(int r = minRmem; r < maxRmem; r++){

          if(!optimization){
            looktableR = r - minRmem - lagRmem;
          } else {
            looktableR = r - minRmem;
          }

          for(int c = minCmem; c < maxCmem; c++){

            if(!optimization){
              looktableC = c - minCmem - lagCmem;
            } else {
              looktableC = c - minCmem;
            }

            weightR = r_memoryKernel.vals[looktableR][looktableC];
            weightW = w_memoryKernel.vals[looktableR][looktableC];

            // reverse "full decay"
            Arena.arrayMemoriesRef[r][c] = Arena.arrayMemoriesRef[r][c] / memoryRD_cplm;
            Arena.arrayMemoriesWork[r][c] = Arena.arrayMemoriesWork[r][c] / memoryWD_cplm;

            Arena.arrayMemoriesRef[r][c] = Arena.arrayMemoriesRef[r][c] -
              (1 - weightR) * Arena.arrayMemoriesRef[r][c] * memoryRD[0] +
              weightR * memoryRL[0];

            Arena.arrayMemoriesWork[r][c] = Arena.arrayMemoriesWork[r][c] -
              (1 - weightW) * Arena.arrayMemoriesWork[r][c] * memoryWD[0] +
              weightW * memoryWL[0];
          }
        }

        // If next time step has valid coordinates
        // seems to be optimization only, check redundancy with
        // kernel check above
        if(optimization){
          if(nextCol != -9999){

            // 2. Calculate movement probability
            sum_weights = 0.0;
            for(int r = minR; r < maxR; r++){
              looktableR = r - minR;

              for(int c = minC; c < maxC; c++){
                looktableC = c - minC;

                familiarity = Arena.arrayMemoriesRef[r][c] - Arena.arrayMemoriesWork[r][c];

                Arena.arrayAttractionWeight[r][c] = stepLengthKernel.vals[looktableR][looktableC] *
                  (Arena.arrayResourceSelection[r][c] * (familiarity + 1));

                sum_weights = sum_weights + Arena.arrayAttractionWeight[r][c];
              }
            }

            // 3. Calculate step likelihood
            if(sum_weights > 0){
              Traj.likelihood[rowCount] = log(Arena.arrayAttractionWeight[nextRow][nextCol] / sum_weights);
            } else {
              Traj.likelihood[rowCount] = log(Arena.arrayAttractionWeight[nextRow][nextCol]);
            }

          }// nextCol section
        } // only during optimization (if clause)
      } // currentCol section
    }  // individual

    // PART 2: SIMULATIONS!
    if(!optimization){
      rowCount = rowCount + 1;

      totalTraj.animalId.push_back(TrajMetrics.animalId[ind]);
      totalTraj.run.push_back(0);
      totalTraj.col.push_back(Traj.col[rowCount]);
      totalTraj.row.push_back(Traj.row[rowCount]);

      totalCount = totalCount + 1;

      // Save info for re-initializing runs
      releaseCol = totalTraj.col[totalCount];
      releaseRow = totalTraj.row[totalCount];

      double** memoriesRefIni;
      double** memoriesWorkIni;
      initialize2D(memoriesRefIni,Arena.nRows,Arena.nCols);
      initialize2D(memoriesWorkIni,Arena.nRows,Arena.nCols);

      for(int r=0;r<Arena.nRows;r++){
        for(int c=0;c<Arena.nCols;c++){
          memoriesRefIni[r][c] = Arena.arrayMemoriesRef[r][c];
          memoriesWorkIni[r][c] = Arena.arrayMemoriesWork[r][c];
        }
      }

      // LOOPS
      for(int run = 1; run <= nSimulatedRuns; run++ ){
        for(ite = 0; ite < nSimulatedSteps; ite++ ){
          // Ensures that the first point is the release coordinate
          // not the last point of the previous run!
          if(ite==0){

            currentCol = releaseCol;
            currentRow = releaseRow;

            // re-update memory
            for(int r=0;r<Arena.nRows;r++){
              for(int c=0;c<Arena.nCols;c++){
                Arena.arrayMemoriesRef[r][c] = memoriesRefIni[r][c];
                Arena.arrayMemoriesWork[r][c] = memoriesWorkIni[r][c];
              }
            }
          } else {

            currentCol = totalTraj.col[totalCount];
            currentRow = totalTraj.row[totalCount];
          }

          minR = currentRow - stepLengthKernel.nCells;
          maxR = currentRow + stepLengthKernel.nCells + 1;
          minC = currentCol - stepLengthKernel.nCells;
          maxC = currentCol + stepLengthKernel.nCells + 1;
          minRmem = currentRow - r_memoryKernel.nCells;
          maxRmem = currentRow + r_memoryKernel.nCells + 1;
          minCmem = currentCol - r_memoryKernel.nCells;
          maxCmem = currentCol + r_memoryKernel.nCells + 1;

          // Corrections to make sure the kernels do not go beyond the study area
          lagR=0;
          lagC=0;

          if(minR<0){lagR=minR;minR=0;}
          if(maxR>Arena.nRows){maxR=Arena.nRows;}
          if(minC<0){lagC=minC;minC=0;}
          if(maxC>Arena.nCols){maxC=Arena.nCols;}
          lagRmem=0;lagCmem=0;
          if(minRmem<0){lagRmem=minRmem;minRmem=0;}
          if(maxRmem>Arena.nRows){maxRmem=Arena.nRows;}
          if(minCmem<0){lagCmem=minCmem;minCmem=0;}
          if(maxCmem>Arena.nCols){maxCmem=Arena.nCols;}

          // 1. Memory dynamics
          for(int r=0;r<Arena.nRows;r++){
            for(int c=0;c<Arena.nCols;c++){
              Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]*memoryRD_cplm;
              Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]*memoryWD_cplm;
            }
          }

          // At cells within neighborhood...
          for(int r=minRmem;r<maxRmem;r++){
            looktableR=r-minRmem-lagRmem;

            for(int c=minCmem;c<maxCmem;c++){
              looktableC=c-minCmem-lagCmem;

              weightR=r_memoryKernel.vals[looktableR][looktableC];
              weightW=w_memoryKernel.vals[looktableR][looktableC];

              Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]/memoryRD_cplm;
              Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]/memoryWD_cplm;

              Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]-
                (1-weightR)*Arena.arrayMemoriesRef[r][c]*memoryRD[0]+
                weightR*memoryRL[0];

              Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]-
                (1-weightW)*Arena.arrayMemoriesWork[r][c]*memoryWD[0]+
                weightW*memoryWL[0];
            }
          }

          // 2. Calculate movement probability
          sum_weights=0.0;
          for(int r=minR;r<maxR;r++){
            looktableR=r-minR-lagR;

            for(int c=minC;c<maxC;c++){
              looktableC=c-minC-lagC;

              familiarity=Arena.arrayMemoriesRef[r][c]-Arena.arrayMemoriesWork[r][c];
              Arena.arrayAttractionWeight[r][c]=stepLengthKernel.vals[looktableR][looktableC]*
                (Arena.arrayResourceSelection[r][c]*(familiarity+1));

              sum_weights = sum_weights + Arena.arrayAttractionWeight[r][c];
            }
          }

          // only use movement kernel
          if(sum_weights<=0){
            for(int r=minR;r<maxR;r++){
              looktableR=r-minR;
              for(int c=minC;c<maxC;c++){
                looktableC=c-minC;
                Arena.arrayAttractionWeight[r][c]=stepLengthKernel.vals[looktableR][looktableC];
                sum_weights=sum_weights+Arena.arrayAttractionWeight[r][c];
              }
            }
          }

          // 3. Calculate random step
          randNumber=R::runif(0,1) * sum_weights;

          // reset value every time
          sumAttractionVectors = 0;

          // evaluate all positions
          for(int r=minR;r<maxR;r++){
            for(int c=minC;c<maxC;c++){
              if(Arena.arrayAttractionWeight[r][c] > 0){
                if(randNumber > sumAttractionVectors){
                  if(randNumber <= (sumAttractionVectors + Arena.arrayAttractionWeight[r][c])){
                    nextCol = c;
                    nextRow = r;
                  }
                }
                sumAttractionVectors = sumAttractionVectors + Arena.arrayAttractionWeight[r][c];
              }
            }
          }

          // 4. Write simulated point
          totalTraj.animalId.push_back(TrajMetrics.animalId[ind]);
          totalTraj.run.push_back(run);
          totalTraj.col.push_back(nextCol);
          totalTraj.row.push_back(nextRow);

          // keep track of how many points written
          totalCount = totalCount + 1;

        }
      }
    } else {
      rowCount = rowCount + 1;
    }

  } // number of animals loop

  // cleanup data frames
  // this avoids runaway memory use on iterations which is
  // especially important during optimization
  // (e.g. ~5GB and climbing for 300 its while constant at about
  // 1.5GB with cleanup)
  arena_cleanup(Arena);
  lookuptable_cleanup(stepLengthKernel);
  lookuptable_cleanup(r_memoryKernel);
  lookuptable_cleanup(w_memoryKernel);


  if(verbose){
    end = std::chrono::high_resolution_clock::now();
    dur = end - begin;
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();

    time_t endTime = time(0);
    Rcout << "gros calculs, run time: "<< ms << " ms" << std::endl;
    Rcout << "MODEL RUN ENDS, run time: "<< endTime-iniTime << " seconds" << endl;
  }

  // split output
  if(optimization){
    return List::create(
      _["likelihood"] = Traj.likelihood
    );
  } else {
    return List::create(
      _["locations"] = DataFrame::create(
        _["id"] = totalTraj.animalId,
        _["run"] = totalTraj.run,
        _["col"] = totalTraj.col,
        _["row"] = totalTraj.row
      ),
      _["resources"] = ListMatrix::create(res)
    );
  }

}
