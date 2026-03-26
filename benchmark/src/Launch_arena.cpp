/*
 * Launch_arena.cpp
 *
 *  Created on: 3 November 2019
 *      Author: Nathan Ranc
 */

#include "header.h"

/*
 *
 * SUMMARY OF ARRAYS IN C++
 *
 */

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



///////////////////////////////////////////////////////////



// Reads the spatial rasters and creates the spatial environment
ArraysDynamics launchArena(std::string raster_directory, std::string output_directory, std::vector<std::string>variable_names, std::vector<double>selection_coef, bool writing_outputs)
{
    ArraysDynamics returnValues;

    // INITIALIZE VARIABLES

    std::ofstream ArenaLogFile;
    if(writing_outputs==true)
    {
            std::string filePath_Alog = output_directory + "arenaLogFile.txt";
            ArenaLogFile.open(filePath_Alog.c_str());
            ArenaLogFile << "ARENA BEGINS" << std::endl;ArenaLogFile<<"\n\n\n";
    }

    std::string filePath, line, extracted;
    int v,c,r;
    double val;
    std::string::size_type sz;


    // INSERT ARENA PARAMETERS
    filePath = raster_directory +  variable_names[0] + ".asc";

    returnValues.nVariables=variable_names.size();
    std::ifstream myfile (filePath.c_str());
    getline (myfile,line);
    line.erase(0,6);
    returnValues.nCols=std::stoi(line);
    getline (myfile,line);
    line.erase(0,6);
    returnValues.nRows=std::stoi(line);
    getline (myfile,line);
    line.erase(0,10);
    returnValues.minXArena=std::stof(line);
    getline (myfile,line);
    line.erase(0,10);
    returnValues.minYArena=std::stof(line);
    getline (myfile,line);
    line.erase(0,9);
    returnValues.resolution=std::stof(line);
    returnValues.maxXArena=returnValues.resolution* double (returnValues.nCols)+returnValues.minXArena;
    returnValues.maxYArena=returnValues.resolution* double (returnValues.nRows)+returnValues.minYArena;
    getline (myfile,line);
    line.erase(0,13);
    returnValues.noDataVal=std::stof(line);

    if(writing_outputs==true)
    {
        ArenaLogFile<<"The Arena has for nRows (y cells): "<< returnValues.nRows<<"\n";
        ArenaLogFile<<"The Arena has for nCols (x cells): "<< returnValues.nCols<<"\n";
        ArenaLogFile<<"The Arena has for resolution: "<< returnValues.resolution<<"\n";
        ArenaLogFile<<"The Arena has for X-corner: "<< returnValues.minXArena<<"\n";
        ArenaLogFile<<"The Arena has for Y-corner: "<< returnValues.minYArena<<"\n";
        ArenaLogFile<<"The Arena has for maximum X: "<< returnValues.maxXArena<<"\n";
        ArenaLogFile<<"The Arena has for maximum Y: "<< returnValues.maxYArena<<"\n";
        ArenaLogFile<<"The Arena has for noDataValue : "<< returnValues.noDataVal<<"\n";
    }


    ///////////////////////////////////////////////////////////


    // INITIALIZE ARENA STRUCTURE
    auto begin = std::chrono::high_resolution_clock::now();
    initialize2D(returnValues.arrayMemoriesRef,returnValues.nRows,returnValues.nCols);
    initialize2D(returnValues.arrayMemoriesWork,returnValues.nRows,returnValues.nCols);
    initialize2D(returnValues.arrayResourceSelection,returnValues.nRows,returnValues.nCols);
    initialize2D(returnValues.arrayAttractionWeight,returnValues.nRows,returnValues.nCols);
    auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << "launchArena -- initialize 2D, run time: "<< ms << " ms" << std::endl;


    ///////////////////////////////////////////////////////////


    // LOAD RASTERS + GET OVERALL RESOURCE VALUE
    begin = std::chrono::high_resolution_clock::now();
    for(v=0;v<returnValues.nVariables;v++)
    {
        // READ RASTERS
        filePath = raster_directory +  variable_names[v] + ".asc";
        if(writing_outputs==true)
        {
            ArenaLogFile<<"Reading " << filePath<<"\n";
        }
        std::ifstream myfile (filePath.c_str());
        getline (myfile,line);
        line.erase(0,6);
        int nCols=std::stoi(line);
        getline (myfile,line);
        line.erase(0,6);
        int nRows=std::stoi(line);
        getline (myfile,line);
        line.erase(0,10);
        getline (myfile,line);
        line.erase(0,10);
        getline (myfile,line);
        line.erase(0,9);
        getline (myfile,line);
        line.erase(0,13);

        if(v==0)
        {
            for(r=0;r<nRows;r++)
            {
                getline (myfile,line);
                std::stringstream ss(line);
                std::string extracted;

                for(c=0;c<nCols;c++)
                {
                		getline(ss, extracted, ' ');
                		val = std::stof (extracted,&sz);
                    if(val==returnValues.noDataVal)
                    {
                        returnValues.arrayResourceSelection[r][c]=0;
                    }
                    else
                    {
                        returnValues.arrayResourceSelection[r][c]=exp(val*selection_coef[v]);
                    }
                }
            }
        }
        else
        {
            for(r=0;r<nRows;r++)
            {
                getline (myfile,line);
                std::stringstream ss(line);
                std::string extracted;

                for(c=0;c<nCols;c++)
                {
                		getline(ss, extracted, ' ');
                		val = std::stof (extracted,&sz);
                		if(val==returnValues.noDataVal)
                		{
                			returnValues.arrayResourceSelection[r][c]=0;
                		}
                		else
                		{
                    		returnValues.arrayResourceSelection[r][c]*=exp(val*selection_coef[v]);
                		}
                }
            }
        }

        if(writing_outputs==true)
        {
            ArenaLogFile<<"Raster layer added!\n";
        }

    }
    end = std::chrono::high_resolution_clock::now();
    dur = end - begin;
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << "launchArena -- read in variables, run time: "<< ms << " ms" << std::endl;


    // WRITE RASTERS
    if(writing_outputs==true)
    {
        std::string filePath_writing = output_directory + "global_resource.asc";
        std::ofstream outputFile2;
        outputFile2.open(filePath_writing.c_str());
        outputFile2 << "ncols "<< returnValues.nCols << std::endl;
        outputFile2 << "nrows "<<returnValues.nRows<<std::endl;
        outputFile2 << "xllcorner "<<returnValues.minXArena<<std::endl;
        outputFile2 << "yllcorner "<<returnValues.minYArena<<std::endl;
        outputFile2 << "cellsize "<<returnValues.resolution<<std::endl;
        outputFile2 << "NODATA_value "<<returnValues.noDataVal<<std::endl<<std::endl;

        for(r=0;r<returnValues.nRows;r++)
        {
        		for(c=0;c<returnValues.nCols;c++)
        		{
        			outputFile2<<returnValues.arrayResourceSelection[r][c]<<" ";
        		}
        }
        outputFile2 <<""<<std::endl;
        outputFile2.close();

        ArenaLogFile<<"\n\n\n";
        ArenaLogFile << "ARENA SUCCESS" << std::endl;
    }

    return returnValues;
}
