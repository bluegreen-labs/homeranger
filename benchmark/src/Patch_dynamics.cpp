/*
 * Patch_dynamics.cpp
 *
 *  Created on: 3 November 2019
 *      Author: Nathan Ranc
 */

#include "header.h"

// Resets the spatial landscapes (memory and attraction) -- to be ran between animals
void arena_renewal(ArraysDynamics & input_arrays, int focPatchX, int focPatchY)
{
    int i, j;

    input_arrays.randmovement=false;
    input_arrays.movement=true;
    input_arrays.checkedPatchCoordinates[0]=focPatchX;
    input_arrays.checkedPatchCoordinates[1]=focPatchY;

    for(i=0;i<input_arrays.nRows;i++)
    {
        for(j=0;j<input_arrays.nCols;j++)
        {
            input_arrays.arrayMemoriesRef[i][j]=0;
            input_arrays.arrayMemoriesWork[i][j]=0;
            input_arrays.arrayAttractionWeight[i][j]=0.0;
        }
    }
}
