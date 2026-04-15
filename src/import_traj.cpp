#include "header.h"

// Imports the trajectories -- transform continuous coordinates into discrete patch ID
structTrajectory launchTrajectoryCoordinatesMat(
        arma::mat locations,
        double resolution,
        double min_x,
        double min_y,
        int n_row,
        int n_col,
        int n_cells_mem
){
    structTrajectory returnValues;
    int counter=0;
    int reserve_number=0;
    double col,row;
    double ratio_position_x, ratio_position_y;
    int ratio_int_position_x, ratio_int_position_y;
    int ind=0;

    bool new_animal;

    // split out the data
    arma::vec id = locations.col(0);
    arma::vec cols = locations.col(1);
    arma::vec rows = locations.col(2);

    int total_relocations = id.size();

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

    for(int r=0; r < total_relocations; r++){
        counter=counter+1;

        if(counter>0){
            returnValues.likelihood.push_back(-9999);
            if(id[r] != ind){
                new_animal=true;
            }else{
                new_animal=false;
            }

            ind = id[r];
            returnValues.animalId.push_back(ind);

            col=cols[r];
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

            row=rows[r];
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
                returnValues.minColMem[counter-1]=returnValues.col[counter-1]-n_cells_mem;
                returnValues.maxColMem[counter-1]=returnValues.col[counter-1]+n_cells_mem;
                returnValues.minRowMem[counter-1]=returnValues.row[counter-1]-n_cells_mem;
                returnValues.maxRowMem[counter-1]=returnValues.row[counter-1]+n_cells_mem;
            }
            else
            {
                if(col!=-9999)
                {

                    if((returnValues.col[counter-1]-n_cells_mem)<returnValues.minColMem[counter-2])
                    {
                        returnValues.minColMem[counter-1]=returnValues.col[counter-1]-n_cells_mem;
                    }
                    else
                    {
                        returnValues.minColMem[counter-1]=returnValues.minColMem[counter-2];
                    }

                    if((returnValues.col[counter-1]+n_cells_mem)>returnValues.maxColMem[counter-2])
                    {
                        returnValues.maxColMem[counter-1]=returnValues.col[counter-1]+n_cells_mem;
                    }
                    else
                    {
                        returnValues.maxColMem[counter-1]=returnValues.maxColMem[counter-2];
                    }

                    if((returnValues.row[counter-1]-n_cells_mem)<returnValues.minRowMem[counter-2])
                    {
                        returnValues.minRowMem[counter-1]=returnValues.row[counter-1]-n_cells_mem;
                    }
                    else
                    {
                        returnValues.minRowMem[counter-1]=returnValues.minRowMem[counter-2];
                    }

                    if((returnValues.row[counter-1]+n_cells_mem)>returnValues.maxRowMem[counter-2])
                    {
                        returnValues.maxRowMem[counter-1]=returnValues.row[counter-1]+n_cells_mem;
                    }
                    else
                    {
                        returnValues.maxRowMem[counter-1]=returnValues.maxRowMem[counter-2];
                    }
                }
                else
                {
                    returnValues.minColMem[counter-1]=returnValues.minColMem[counter-2];
                    returnValues.maxColMem[counter-1]=returnValues.maxColMem[counter-2];
                    returnValues.minRowMem[counter-1]=returnValues.minRowMem[counter-2];
                    returnValues.maxRowMem[counter-1]=returnValues.maxRowMem[counter-2];
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
