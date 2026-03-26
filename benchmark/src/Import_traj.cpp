/*
 * Import_traj.cpp
 *
 *  Created on: 3 November 2019
 *      Author: Nathan Ranc
 */

#include "header.h"

// Imports the trajectories -- transform continuous coordinates into discrete patch ID
structTrajectory launchTrajectoryCoordinates(std::string path, double resolution, double min_x, double min_y, int n_row, int n_col, int n_cells_mem)
{
    structTrajectory returnValues;
    int counter=0;
    int total_relocations=0;
    int reserve_number=0;
    double col,row;
    double ratio_position_x, ratio_position_y;
    int ratio_int_position_x, ratio_int_position_y;
    int ind=0;

    bool new_animal;

    std::string line, extracted;
    std::string filePath_sLength = path;
    std::ifstream sLength (filePath_sLength.c_str());

    // Count number of relocations in file
    while(getline(sLength, line))
    {
        total_relocations=total_relocations+1;
    }
    total_relocations=total_relocations-1;

    // Initialize
    reserve_number=total_relocations+5;
    returnValues.animalId.reserve(reserve_number);
    returnValues.col.reserve(reserve_number);
    returnValues.row.reserve(reserve_number);
    returnValues.minColMem.reserve(reserve_number);
    returnValues.maxColMem.reserve(reserve_number);
    returnValues.minRowMem.reserve(reserve_number);
    returnValues.maxRowMem.reserve(reserve_number);
    returnValues.resourceSelection.reserve(reserve_number);
    returnValues.likelihood.reserve(reserve_number);

    // Read file
    std::string line2, extracted2;
    std::string::size_type sz2;      // alias of size_t
    std::ifstream sLength2 (filePath_sLength.c_str());

    while(getline(sLength2, line2))
    {
        counter=counter+1;

        if(counter>1)
        {
            returnValues.likelihood.push_back(-9999);

            extracted2 = line2.substr(0, line2.find(","));
            if((std::stoi (extracted2,&sz2))!=ind){new_animal=true;}
            else{new_animal=false;}
            ind=std::stoi (extracted2,&sz2);
            returnValues.animalId.push_back(ind);
            line2=line2.erase(0, line2.find(",")+1);

            extracted2 = line2.substr(0, line2.find(","));
            col=std::stod (extracted2,&sz2);
            if(col==-9999)
            {
                returnValues.col.push_back((int)col);
            }
            else
            {
                ratio_position_x=(col-min_x)/resolution;
                ratio_int_position_x=(int) ratio_position_x;
                returnValues.col.push_back(ratio_int_position_x);
            }
            line2=line2.erase(0, line2.find(",")+1);

            extracted2 = line2.substr(0, line2.find(" "));
            row=std::stod (extracted2,&sz2);
            if(row==-9999)
            {
                returnValues.row.push_back((int)row);
            }
            else
            {
                ratio_position_y=(row-min_y)/resolution;;
                ratio_int_position_y=(int) ceil(ratio_position_y);
                returnValues.row.push_back(n_row-ratio_int_position_y);
            }

            // Definition of the animal's bounding box
            if(new_animal==true)
            {
            		returnValues.minColMem[counter-2]=returnValues.col[counter-2]-n_cells_mem;
            		returnValues.maxColMem[counter-2]=returnValues.col[counter-2]+n_cells_mem;
            		returnValues.minRowMem[counter-2]=returnValues.row[counter-2]-n_cells_mem;
            		returnValues.maxRowMem[counter-2]=returnValues.row[counter-2]+n_cells_mem;
            }
            else
            {
            		if(col!=-9999)
            		{

                		if((returnValues.col[counter-2]-n_cells_mem)<returnValues.minColMem[counter-3])
                		{
                    		returnValues.minColMem[counter-2]=returnValues.col[counter-2]-n_cells_mem;
                		}
                		else
                		{
                    		returnValues.minColMem[counter-2]=returnValues.minColMem[counter-3];
                		}

                		if((returnValues.col[counter-2]+n_cells_mem)>returnValues.maxColMem[counter-3])
                		{
                    		returnValues.maxColMem[counter-2]=returnValues.col[counter-2]+n_cells_mem;
                		}
                		else
                		{
                    		returnValues.maxColMem[counter-2]=returnValues.maxColMem[counter-3];
                		}

                		if((returnValues.row[counter-2]-n_cells_mem)<returnValues.minRowMem[counter-3])
                		{
                    		returnValues.minRowMem[counter-2]=returnValues.row[counter-2]-n_cells_mem;
                		}
                		else
                		{
                    		returnValues.minRowMem[counter-2]=returnValues.minRowMem[counter-3];
                		}

                		if((returnValues.row[counter-2]+n_cells_mem)>returnValues.maxRowMem[counter-3])
                		{
                    		returnValues.maxRowMem[counter-2]=returnValues.row[counter-2]+n_cells_mem;
                		}
                		else
                		{
                    		returnValues.maxRowMem[counter-2]=returnValues.maxRowMem[counter-3];
                		}
            		}
            		else
            		{
                		returnValues.minColMem[counter-2]=returnValues.minColMem[counter-3];
                		returnValues.maxColMem[counter-2]=returnValues.maxColMem[counter-3];
                		returnValues.minRowMem[counter-2]=returnValues.minRowMem[counter-3];
                		returnValues.maxRowMem[counter-2]=returnValues.maxRowMem[counter-3];
            		}
            }
        }
    }

    return returnValues;
}



// Gets the trajectory summary metrics
structSummaryTraj getTrajectoryMetrics(structTrajectory traj)
{
    int l;
    int ind=-9999;
    int indCount=0;
    int count;
    structSummaryTraj returnValues;

    returnValues.totalLength=traj.animalId.size();

    returnValues.animalId.reserve(100);
    returnValues.individualCount.reserve(100);

    for(l=0;l<returnValues.totalLength;l++)
    {
        if(traj.animalId[l]!=ind)
        {
            indCount=indCount+1;
            ind=traj.animalId[l];

            returnValues.animalId.push_back(traj.animalId[l]);

            if(indCount!=1)
            {
                returnValues.individualCount.push_back(count);
            }

            count=1;
        }

        else
        {
            count=count+1;
        }
    }

    returnValues.individualCount.push_back(count);

    return returnValues;
}





// Gets the trajectory's associated spatial information
    // Extracts covariate values for all fixes
    // Note.: this might give an error if the points are outside the rasters!
void joinTrajectoryLandscape(std::string output_directory, structTrajectory & loadedTraj, double** resource_selection_array)
{
    //time_t iniTime = time(0);

    int nRelocations=loadedTraj.col.size();

    for(int r=0;r<nRelocations;r++)
    {
        if(loadedTraj.col[r]==-9999)
        {
            loadedTraj.resourceSelection[r]=0;
        }
        else
        {
            loadedTraj.resourceSelection[r]=resource_selection_array[loadedTraj.row[r]][loadedTraj.col[r]];
        }
    }

    std::ofstream outputFile2;

    std::string filePath_writing = output_directory + "trajectory_covariates.csv";
    outputFile2.open(filePath_writing.c_str());
    outputFile2 << std::fixed;
    outputFile2 << "animal_id,r_patch,c_patch,resource_selection"<<std::endl;
    for(int r=0;r<nRelocations;r++)
    {
        outputFile2 << loadedTraj.animalId[r] <<","<<
        			loadedTraj.row[r] <<","<<loadedTraj.col[r]<<","<<
                loadedTraj.resourceSelection[r]<<std::endl;                                                                                                      ;
    }

    outputFile2.close();
}

