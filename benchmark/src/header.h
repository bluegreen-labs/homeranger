/*
 * header.h
 *
 *  Created on: 3 November 2019
 *      Author: Nathan Ranc
 */

#ifndef HEADER_H_
#define HEADER_H_



/*
 *
 * LIBRARIES
 *
 */

#include    <stdlib.h>
#include    <cmath>
#include    <math.h>
#include    <iostream>
#include    <fstream>
#include    <string>
#include    <random>
#include    <sys/types.h>
#include    <dirent.h>
#include    <errno.h>
#include    <vector>
#include    <ctime>
#include    <chrono>
#include    <sstream>
#include    <algorithm>


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


/*
 *
 * STRUCTURES
 *
 */

struct structConfig
{
    std::string outputDirectoryPath;             // output directory path
    std::string rasterDirectoryPath;             // raster directory path
    std::string trajectoryPath;                  // movement trajectories (regularized, with interpolated NAs)

    bool writingOututs=false;                    // write test rasters

    std::string modelType;                       // "kernel" or "simulation"

    double thresholdApproxKernel;                // distance used to threshold attraction calculation
    double thresholdMemoryKernel;                // distance used to threshold attraction calculation

    std::vector<std::string>resourceNames;       // list of resource names
    std::vector<double>selectionCoef;            // corresponding variable selection coefficient

    double memoryRL;                             // learning rate for reference memory
    double memoryWL;                             // learning rate for working memory  >0 (because division by this 1-parameter in code)
    double memoryRD;                             // decay rate for reference memory	>0 (because division by this 1-parameter in code)
    double memoryWD;                             // decay rate for working memory
    double memoryRD_cplm;                        // decay rate for reference memory -- complement (1-d)
    double memoryWD_cplm;                        // decay rate for working memory -- complement (1-d)

    double memoryRDist;                          // R mem exponential decay rate for distance devaluation
    double memoryWDist;                          // W mem exponential decay rate for distance devaluation >0
    double stepLengthDist;                       // step length exponential decay rate for distance devaluation
    double stepLengthShape;						// shape parameter for

    int nSimulatedSteps;							// number of simulated steps following observed locations
    int nSimulatedRuns;							// number of simulated runs
};


struct ArraysDynamics
{
    double** arrayMemoriesRef;               // raster of reference memory of any given cell
    double** arrayMemoriesWork;              // raster of working memory of any given cell
    double** arrayResourceSelection;         // raster of combined resource selection of any given cell
    double** arrayAttractionWeight;          // raster of attraction weight of any given cell (normalized)

    int nRows;                              	// number of rows on an .asc file
    int nCols;                              	// number of columns on an .asc file
    double resolution;
    double minXArena;
    double maxXArena;
    double minYArena;
    double maxYArena;
    double noDataVal;
    int nVariables;

    int usedPatchCoordinates[2];
    int checkedPatchCoordinates[2];

    bool movement;
    bool randmovement;
};


struct lookupTable
{
    int nCells;
    double** vals;
};


struct structTrajectory
{
    std::vector<int>animalId;
    std::vector<int>col;
    std::vector<int>row;
    std::vector<int>minColMem;
    std::vector<int>maxColMem;
    std::vector<int>minRowMem;
    std::vector<int>maxRowMem;
    std::vector<double>resourceSelection;
    std::vector<double>likelihood;
};

struct structTrajectorySimul
{
    std::vector<int>animalId;
    std::vector<int>run;
    std::vector<int>col;
    std::vector<int>row;
    std::vector<int>minColMem;
    std::vector<int>maxColMem;
    std::vector<int>minRowMem;
    std::vector<int>maxRowMem;
};


struct structSummaryTraj
{
    int totalLength;
    std::vector<int>animalId;
    std::vector<int>individualCount;
};


struct position
{
    double x;
    double y;
};


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


/*
 *
 * FUNCTIONS
 *
 */

// Configuration.cpp
structConfig launchConfig_txt_sep(std::string path_to_config);

// Writing_outputs.cpp
void writeTrajCompoToLog(std::ofstream & logfile, structTrajectory traj, structSummaryTraj traj_metrics);

// Launch_arena.cpp
void initialize2D(double** & vector, int nR, int nC);
ArraysDynamics launchArena(std::string raster_directory, std::string output_directory, std::vector<std::string>variable_names, std::vector<double>selection_coef, bool writing_outputs);
void initialize2D_call(double** & vector, int nR, int nC);

// Import_traj.cpp
structTrajectory launchTrajectoryCoordinates(std::string path, double resolution, double min_x, double min_y, int n_row, int n_col, int n_cells_mem);
structSummaryTraj getTrajectoryMetrics(structTrajectory loadedTraj);
void joinTrajectoryLandscape(std::string output_directory, structTrajectory & loadedTraj, double** resource_selection_array);

// Dist_lookup_table.cpp
lookupTable iniApproxKernel (double distance_threshold, double resolution, double spatial_decay);
lookupTable iniApproxKernelStepLength (double distance_threshold, double resolution, double spatial_decay, double shape, double residence_p);

// Patch_dynamics.cpp
void arena_renewal(ArraysDynamics & inputArrays, int focPatchX, int focPatchY);

// Likelihood.cpp
void writeObjFunction(structTrajectory traj, std::string output_directory, bool write_outputs);


#endif /* HEADER_H_ */
