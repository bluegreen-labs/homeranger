/*
 * Likelihood.cpp
 *
 *  Created on: 3 November 2019
 *      Author: Nathan Ranc
 */

#include "header.h"

// Returns the objective function value
void writeObjFunction(structTrajectory traj, std::string output_directory, bool write_outputs)
{
    double objectiveFunction=0.0;

    for(int i=0;i<traj.likelihood.size();i++)
    {
            if(traj.likelihood[i]!=-9999)
            {
                objectiveFunction=objectiveFunction+traj.likelihood[i];
            }
    }

    std::ofstream objFunctionFile;
    std::string filePath_output =  output_directory + "objective_function.csv";
    objFunctionFile.setf(std::ios::fixed, std::ios::floatfield);
    objFunctionFile.precision(7);
    objFunctionFile.open(filePath_output.c_str());
    objFunctionFile << objectiveFunction<< std::endl;
    objFunctionFile.close();

    if(write_outputs==true)
    {
        std::ofstream objFunctionDetailFile;
        std::string filePath_output2 =  output_directory + "objective_function_detail.csv";
        objFunctionDetailFile.setf(std::ios::fixed, std::ios::floatfield);
        objFunctionDetailFile.precision(7);
        objFunctionDetailFile.open(filePath_output2.c_str());
        objFunctionDetailFile << "animal_id,r_patch,c_patch,likelihood"<< std::endl;

        for(int i=0;i<traj.likelihood.size();i++)
        {
            objFunctionDetailFile<<traj.animalId[i]<<","<<traj.row[i]<<","<<traj.col[i]<<","<<traj.likelihood[i]<< std::endl;
        }

        objFunctionDetailFile.close();
    }
}
