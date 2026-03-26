/*
 * Writing_outputs.cpp
 *
 *  Created on: 3 November 2019
 *      Author: Nathan Ranc
 */

#include "header.h"

// Opens file
void openFile(std::ofstream & file, std::string fileName, std::string PrRrep, std::string PrWrep, std::string dir)
{
    std::string filePath = dir + "outputs/" + fileName + "_refMem" + PrRrep + "_workMem" + PrWrep + ".csv";
    file.open(filePath.c_str());
}

// Closes file
void closeFile(std::ofstream & file)
{
    file.close();
}

// Saves trajectory summary
void writeTrajCompoToLog(std::ofstream & logfile, structTrajectory traj, structSummaryTraj traj_metrics)
{
    int i;
    int nAnimals=traj_metrics.animalId.size();

    logfile << "TRAJECTORY LOADED" << std::endl;
    logfile << "Number of animals: "<< nAnimals << std::endl;
    logfile << "Number of relocations: "<< traj.animalId.size() << std::endl;

    for(i=0;i<nAnimals;i++)
    {
        logfile << "Animals ID: "<< traj_metrics.animalId[i] << ", N locs.: "<< traj_metrics.individualCount[i] << std::endl;
    }
    logfile << "" << std::endl;
}
