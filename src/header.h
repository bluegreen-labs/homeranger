/*
 *  Created on: 3 November 2019
 *      Author: Nathan Ranc
 */

#ifndef HEADER_H_
#define HEADER_H_

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

// structures

struct ArraysDynamics
{
    double** arrayMemoriesRef;
    double** arrayMemoriesWork;
    double** arrayResourceSelection;
    double** arrayAttractionWeight;

    int nRows;
    int nCols;

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


// functions

// arena_dynamics.cpp
void arena_renewal(ArraysDynamics & inputArrays, int focPatchX, int focPatchY);
void arena_cleanup(ArraysDynamics & inputArrays);
void lookuptable_cleanup(lookupTable & inputArrays);

void initialize2D(double** & vector, int nR, int nC);
ArraysDynamics launchArena(
        int nRow,
        int nCol
);
void initialize2D_call(double** & vector, int nR, int nC);

// Import_traj.cpp
structTrajectory launchTrajectoryCoordinates(
        std::string path,
        double resolution,
        double min_x,
        double min_y,
        int n_row,
        int n_col,
        int n_cells_mem
);

structSummaryTraj getTrajectoryMetrics(structTrajectory loadedTraj);

// kernels.cpp
lookupTable iniApproxKernel (
        double distance_threshold,
        double resolution,
        double spatial_decay
);

lookupTable iniApproxKernelStepLength (
        double distance_threshold,
        double resolution,
        double spatial_decay,
        double shape,
        double residence_p
);


#endif /* HEADER_H_ */
