#include "header.h"

// Initialization
void initialize2D(double** & vector, int nR, int nC)
{
  vector = new double*[nR];
  for (int i = 0; i < nR; i++)
  {
    vector[i] = new double[nC];
  }
}

void initialize2D_call(double** & vector, int nR, int nC)
{
  vector = new double*[nR];
  for (int i = 0; i < nR; i++)
  {
    vector[i] = new double[nC];
  }
}

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

// Resets the spatial landscapes (memory and attraction)
void lookuptable_cleanup(lookupTable & dynamics)
{
  // 1. Deallocate the inner arrays (the columns) for each row
  for (int i = 0; i < dynamics.nCells; ++i) {
    delete[] dynamics.vals[i];
  }

  // 2. Deallocate the array of row pointers
  delete[] dynamics.vals;
}

// Resets the spatial landscapes (memory and attraction)
void arena_cleanup(ArraysDynamics & dynamics)
{
  // 1. Deallocate the inner arrays (the columns) for each row
  for (int i = 0; i < dynamics.nRows; ++i) {
    delete[] dynamics.arrayMemoriesWork[i];
    delete[] dynamics.arrayMemoriesRef[i];
    delete[] dynamics.arrayResourceSelection[i];
    delete[] dynamics.arrayAttractionWeight[i];
  }

  // 2. Deallocate the array of row pointers
  delete[] dynamics.arrayMemoriesWork;
  delete[] dynamics.arrayMemoriesRef;
  delete[] dynamics.arrayResourceSelection;
  delete[] dynamics.arrayAttractionWeight;
}

// If the structure itself was allocated on the heap:
// delete dynamicsPtr; // Deallocate the structure itself

// Reads the spatial rasters and creates the spatial environment
ArraysDynamics launchArena(
    int nRows,
    int nCols
){
  ArraysDynamics returnValues;

  returnValues.nRows = nRows;
  returnValues.nCols = nCols;

  initialize2D(returnValues.arrayMemoriesRef,nRows,nCols);
  initialize2D(returnValues.arrayMemoriesWork,nRows,nCols);
  initialize2D(returnValues.arrayResourceSelection,nRows,nCols);
  initialize2D(returnValues.arrayAttractionWeight,nRows,nCols);
  return returnValues;
}
