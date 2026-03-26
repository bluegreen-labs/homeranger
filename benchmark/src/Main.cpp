/*
 * Main.cpp
 *
 *  Created on: 3 November 2019
 *      Author: Nathan Ranc
 */



#include "header.h"

using namespace std;

int main(int argc, char* argv[])
{

    time_t iniTime = time(0);

    bool execute_command_line=true;        	// if execute_command_line=TRUE  => config file via command line
                                        		// if execute_command_line=FALSE => default config file (DEBUGGING)


    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////


    /*
     *
     * CONFIG FILE
     *
     */

    string configFile;

    if(execute_command_line==false)		// obtains default config file
    {
		// RESOURCE-ONLY MODEL
		//configFile="/Path/to/config_best_Mres_fitting.txt";
		//configFile="/Path/to/config_best_Mres_simulation_1y.txt";
		//configFile="/Path/to/config_best_Mres_simulation_2y.txt";

		// MEORY-BASED MODEL
    		//configFile="/Path/to/config_best_Mmem_fitting.txt";
    		//configFile="/Path/to/config_best_Mmem_simulation_1y.txt";
    		//configFile="/Path/to/config_best_Mmem_simulation_2y.txt";
    }
    else		// obtains config file path via command line arguments
    {
        if (argc != 3)
        {
            std::cout << "PROGRAM ENDED: expecting a single parameter (too few or too many)"<<endl;
            exit(0);
        }
        else
        {
            for (int a = 1; a < argc; a++)
            {
                std::string arg = argv[a];
                if (a + 1 != argc)
                {
                    if (arg == "-config")
                    {
                        configFile=argv[a+1];
                    }
                 }
            }
        }
    }

    structConfig configArg=launchConfig_txt_sep(configFile);


    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////


    /*
     *
     * LOG FILE
     *
     */

    std::ofstream MainLogFile;
    if(configArg.writingOututs==true)
    {
        std::string filePath_log = configArg.outputDirectoryPath + "mainLogFile.txt";
        MainLogFile.open(filePath_log.c_str());
        MainLogFile << "MODEL FITTING BEGINS" << endl << endl;
    }


    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////


    /*
     *
     * VARIABLE DECLARATION
     *
     */

    double familiarity, sum_weights;
    int ind, ite, rowCount;
    int currentCol,currentRow,nextCol, nextRow, releaseCol, releaseRow;
    int minR, maxR, minC, maxC;
    int minRmem, maxRmem, minCmem, maxCmem;
    int lagR, lagC, lagRmem, lagCmem;
    int looktableR, looktableC;
    double weightR, weightW;
    int focusPatchX=0, focusPatchY=0;

    if(configArg.writingOututs==true)
    {
        MainLogFile << "Variables declared" << endl;
    }


    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////


    /*
    *
    * ARENA INITIALIZATION
    *
    */

    ArraysDynamics Arena = launchArena(configArg.rasterDirectoryPath,configArg.outputDirectoryPath,configArg.resourceNames,configArg.selectionCoef,configArg.writingOututs);
    if(configArg.writingOututs==true)
    {
            MainLogFile << "Arena initialized" << endl;
    }


    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////


    /*
     *
     * DISTANCE LOOKUP TABLES
     *
     */

    auto begin = std::chrono::high_resolution_clock::now();
    lookupTable stepLengthKernel=iniApproxKernelStepLength(configArg.thresholdApproxKernel, Arena.resolution, configArg.stepLengthDist, configArg.stepLengthShape,0);
    lookupTable r_memoryKernel=iniApproxKernel(configArg.thresholdMemoryKernel, Arena.resolution, configArg.memoryRDist);
    lookupTable w_memoryKernel=iniApproxKernel(configArg.thresholdMemoryKernel, Arena.resolution, configArg.memoryWDist);
    auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << "lookup tables, run time: "<< ms << " ms" << std::endl;


    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////


    /*
     *
     * LOADING TRAJECTORY
     *
     */

    structTrajectory Traj=launchTrajectoryCoordinates(configArg.trajectoryPath, Arena.resolution, Arena.minXArena, Arena.minYArena, Arena.nRows, Arena.nCols, r_memoryKernel.nCells);
    structSummaryTraj TrajMetrics=getTrajectoryMetrics(Traj);
    int nAnimals=TrajMetrics.animalId.size();

    if(configArg.writingOututs==true)
    {
    		writeTrajCompoToLog(MainLogFile, Traj, TrajMetrics);
    		joinTrajectoryLandscape(configArg.outputDirectoryPath, Traj, Arena.arrayResourceSelection);
    }



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
/////////////////////KERNEL FITTING////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    if(configArg.modelType=="kernel")
    {
	    begin = std::chrono::high_resolution_clock::now();
        rowCount=-1;
		lagR=0;lagC=0;lagRmem=0;lagCmem=0;

        // START OF INDIVIDUAL LOOP
        for(ind=0;ind<nAnimals;ind++)
        {
            // Re-initialization of the arena
            arena_renewal(Arena, focusPatchX, focusPatchY);

            // START OF RELOCATION ITERATIONS
            for(ite=0;ite<TrajMetrics.individualCount[ind]-1;ite++)         // RELOCATION LOOP (from release to 2nd from last point)
            {
            		rowCount=rowCount+1;

            		currentCol=Traj.col[rowCount];
            		currentRow=Traj.row[rowCount];
            		nextCol=Traj.col[rowCount+1];
            		nextRow=Traj.row[rowCount+1];

            		minR=currentRow-stepLengthKernel.nCells;
            		maxR=currentRow+stepLengthKernel.nCells+1;
            		minC=currentCol-stepLengthKernel.nCells;
            		maxC=currentCol+stepLengthKernel.nCells+1;
            		minRmem=currentRow-r_memoryKernel.nCells;
            		maxRmem=currentRow+r_memoryKernel.nCells+1;
            		minCmem=currentCol-r_memoryKernel.nCells;
            		maxCmem=currentCol+r_memoryKernel.nCells+1;


                // 1. Memory dynamics
                for(int r=Traj.minRowMem[rowCount];r<=Traj.maxRowMem[rowCount];r++)  // At all cells within animal's bounding box...
                {
                		for(int c=Traj.minColMem[rowCount];c<=Traj.maxColMem[rowCount];c++)
                    {
                			Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]*configArg.memoryRD_cplm;  //"full" memory decay
                			Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]*configArg.memoryWD_cplm;
                    }
                }

                if(currentCol!=-9999)      // If current timestep has valid coordinates
                {
                    for(int r=minRmem;r<maxRmem;r++)      // At cells within neighborhood...
                    {
                        looktableR=r-minRmem;

                        for(int c=minCmem;c<maxCmem;c++)
                        {
                        		looktableC=c-minCmem;
                        		weightR=r_memoryKernel.vals[looktableR][looktableC];
                        		weightW=w_memoryKernel.vals[looktableR][looktableC];

                        		Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]/configArg.memoryRD_cplm;  	// reverse "full decay"
                        		Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]/configArg.memoryWD_cplm;	// reverse "full decay"

                        		Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]-
                                        (1-weightR)*Arena.arrayMemoriesRef[r][c]*configArg.memoryRD+
                                        weightR*configArg.memoryRL;

                        		Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]-
                                        (1-weightW)*Arena.arrayMemoriesWork[r][c]*configArg.memoryWD+
                                        weightW*configArg.memoryWL;
                        }
                    }

                    if(nextCol!=-9999) // If next timestep has valid coordinates
                    {
                        // 2. Calculate movement probability
                        sum_weights=0.0;
                        for(int r=minR;r<maxR;r++)
                        {
                        		looktableR=r-minR;

                        		for(int c=minC;c<maxC;c++)
                            {
                            		looktableC=c-minC;

                        			familiarity=Arena.arrayMemoriesRef[r][c]-Arena.arrayMemoriesWork[r][c];

                                  Arena.arrayAttractionWeight[r][c]=stepLengthKernel.vals[looktableR][looktableC]*
                                          (Arena.arrayResourceSelection[r][c]*(familiarity+1));

                                  sum_weights=sum_weights+Arena.arrayAttractionWeight[r][c];
                            }
                        }


                        // 3. Calculate step likelihood
                        if(sum_weights>0)
                        {
                        		Traj.likelihood[rowCount]=log(Arena.arrayAttractionWeight[nextRow][nextCol]/sum_weights);
                        }
                        else
                        {
                    			Traj.likelihood[rowCount]=log(Arena.arrayAttractionWeight[nextRow][nextCol]);
                        }
                    }
                }
            }

            rowCount=rowCount+1;
        }

        end = std::chrono::high_resolution_clock::now();
        dur = end - begin;
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
        std::cout << "gros calculs, run time: "<< ms << " ms" << std::endl;

        writeObjFunction(Traj,configArg.outputDirectoryPath,configArg.writingOututs);
    }




///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
//////////////////////SIMULATIONS//////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


    if(configArg.modelType=="simulation")
    {
        structTrajectorySimul totalTraj;
        totalTraj.animalId.reserve(TrajMetrics.totalLength+nAnimals*configArg.nSimulatedSteps*configArg.nSimulatedRuns);
        totalTraj.run.reserve(TrajMetrics.totalLength+nAnimals*configArg.nSimulatedSteps*configArg.nSimulatedRuns);
        totalTraj.col.reserve(TrajMetrics.totalLength+nAnimals*configArg.nSimulatedSteps*configArg.nSimulatedRuns);
        totalTraj.row.reserve(TrajMetrics.totalLength+nAnimals*configArg.nSimulatedSteps*configArg.nSimulatedRuns);

        double sumAttractionVectors,randNumber;
		lagR=0;lagC=0;lagRmem=0;lagCmem=0;
        int totalCount;
        rowCount=-1;
        totalCount=-1;

        // START OF INDIVIDUAL LOOP
        for(ind=0;ind<nAnimals;ind++)
        {
            // Re-initialization of the arena
            arena_renewal(Arena, focusPatchX, focusPatchY);


    			// PART 1: UPDATE PRIOR MEMORY!
            if(TrajMetrics.individualCount[ind]>1)
            {
                for(ite=0;ite<TrajMetrics.individualCount[ind]-1;ite++)         // START OF RELOCATION ITERATIONS
                {
                    rowCount=rowCount+1;
                    totalCount=totalCount+1;
                    currentCol=Traj.col[rowCount];
                    currentRow=Traj.row[rowCount];
                    minRmem=currentRow-r_memoryKernel.nCells;
                    maxRmem=currentRow+r_memoryKernel.nCells+1;
                    minCmem=currentCol-r_memoryKernel.nCells;
                    maxCmem=currentCol+r_memoryKernel.nCells+1;

            			// Corrections to make sure the kernels do not go beyond the study area
            			lagRmem=0;lagCmem=0;
            			if(minRmem<0){lagRmem=minRmem;minRmem=0;}
            			if(maxRmem>Arena.nRows){maxRmem=Arena.nRows;}
            			if(minCmem<0){lagCmem=minCmem;minCmem=0;}
            			if(maxCmem>Arena.nCols){maxCmem=Arena.nCols;}


                    // 1. Memory dynamics
            			for(int r=Traj.minRowMem[rowCount];r<=Traj.maxRowMem[rowCount];r++)  // At all cells within animal's bounding box...
            			{
            				for(int c=Traj.minColMem[rowCount];c<=Traj.maxColMem[rowCount];c++)
            				{
            					Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]*configArg.memoryRD_cplm;  //"full" memory decay
            					Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]*configArg.memoryWD_cplm;
                        }
                    }

                    if(currentCol!=-9999)
                    {
                        for(int r=minRmem;r<maxRmem;r++)      // At cells within neighborhood...
                        {
                            looktableR=r-minRmem-lagRmem;
                            // looktableR=r-minRmem;

                            for(int c=minCmem;c<maxCmem;c++)
                            {
                        			looktableC=c-minCmem-lagCmem;
                        			// looktableC=c-minCmem;

                            		weightR=r_memoryKernel.vals[looktableR][looktableC];
                            		weightW=w_memoryKernel.vals[looktableR][looktableC];

                            		Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]/configArg.memoryRD_cplm;  	// reverse "full decay"
                            		Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]/configArg.memoryWD_cplm;	// reverse "full decay"

                            		Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]-
                                            (1-weightR)*Arena.arrayMemoriesRef[r][c]*configArg.memoryRD+
                                            weightR*configArg.memoryRL;

                            		Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]-
                                            (1-weightW)*Arena.arrayMemoriesWork[r][c]*configArg.memoryWD+
                                            weightW*configArg.memoryWL;
                            }
                        }
                    }

                    // 3. Update total trajectory
                    totalTraj.animalId.push_back(TrajMetrics.animalId[ind]);
                    totalTraj.run.push_back(0);
                    totalTraj.col.push_back(currentCol);
                    totalTraj.row.push_back(currentRow);
                }
            }


            // PART 2: SIMULATIONS!
            rowCount=rowCount+1;
            totalCount=totalCount+1;
            totalTraj.animalId.push_back(TrajMetrics.animalId[ind]);
            totalTraj.run.push_back(0);
            totalTraj.col.push_back(Traj.col[rowCount]);
            totalTraj.row.push_back(Traj.row[rowCount]);

            // Save info for re-initilizing runs
            releaseCol=totalTraj.col[totalCount];
			releaseRow=totalTraj.row[totalCount];
		    double** memoriesRefIni;
		    double** memoriesWorkIni;
		    initialize2D_call(memoriesRefIni,Arena.nRows,Arena.nCols);
		    initialize2D_call(memoriesWorkIni,Arena.nRows,Arena.nCols);


            for(int r=0;r<Arena.nRows;r++)
            {
                for(int c=0;c<Arena.nCols;c++)
                {
                		memoriesRefIni[r][c]=Arena.arrayMemoriesRef[r][c];
                		memoriesWorkIni[r][c]=Arena.arrayMemoriesWork[r][c];
                }
            }

            // LOOPS
            for(int run=1;run<=configArg.nSimulatedRuns;run++)		// RUN LOOP
            {
                for(ite=0;ite<configArg.nSimulatedSteps;ite++)		// ITERATION LOOP
                {
                		if(ite==0)	// Ensures that the first point is the release coordinate, not the last point of the previous run!
                		{
                			currentCol=releaseCol;
                			currentRow=releaseRow;

                			// re-update memory
                			for(int r=0;r<Arena.nRows;r++)
                			{
                				for(int c=0;c<Arena.nCols;c++)
                				{
                					Arena.arrayMemoriesRef[r][c]=memoriesRefIni[r][c];
                					Arena.arrayMemoriesWork[r][c]=memoriesWorkIni[r][c];
                				}
                			}
                		}
                		else
                		{
                			currentCol=totalTraj.col[totalCount];
                			currentRow=totalTraj.row[totalCount];
                		}

            			minR=currentRow-stepLengthKernel.nCells;
            			maxR=currentRow+stepLengthKernel.nCells+1;
            			minC=currentCol-stepLengthKernel.nCells;
            			maxC=currentCol+stepLengthKernel.nCells+1;
            			minRmem=currentRow-r_memoryKernel.nCells;
            			maxRmem=currentRow+r_memoryKernel.nCells+1;
            			minCmem=currentCol-r_memoryKernel.nCells;
            			maxCmem=currentCol+r_memoryKernel.nCells+1;

            			// Corrections to make sure the kernels do not go beyond the study area
            			lagR=0;lagC=0;
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
                    for(int r=0;r<Arena.nRows;r++)  			// At all cells...
                    {
                    		for(int c=0;c<Arena.nCols;c++)
                        {
                    			Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]*configArg.memoryRD_cplm;  //"full" memory decay
                    			Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]*configArg.memoryWD_cplm;
                        }
                    }

                    for(int r=minRmem;r<maxRmem;r++)      // At cells within neighborhood...
                    {
                        looktableR=r-minRmem-lagRmem;

                        for(int c=minCmem;c<maxCmem;c++)
                        {
                   			looktableC=c-minCmem-lagCmem;

                   			weightR=r_memoryKernel.vals[looktableR][looktableC];
                        		weightW=w_memoryKernel.vals[looktableR][looktableC];

                        		Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]/configArg.memoryRD_cplm;  	// reverse "full decay"
                        		Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]/configArg.memoryWD_cplm;	// reverse "full decay"

                        		Arena.arrayMemoriesRef[r][c]=Arena.arrayMemoriesRef[r][c]-
                                        (1-weightR)*Arena.arrayMemoriesRef[r][c]*configArg.memoryRD+
                                        weightR*configArg.memoryRL;

                        		Arena.arrayMemoriesWork[r][c]=Arena.arrayMemoriesWork[r][c]-
                                        (1-weightW)*Arena.arrayMemoriesWork[r][c]*configArg.memoryWD+
                                        weightW*configArg.memoryWL;
                        }
                    }

                    // 2. Calculate movement probability
                    sum_weights=0.0;
                    for(int r=minR;r<maxR;r++)
                    {
                   		looktableR=r-minR-lagR;

                    		for(int c=minC;c<maxC;c++)
                    		{
                       		looktableC=c-minC-lagC;

                       		familiarity=Arena.arrayMemoriesRef[r][c]-Arena.arrayMemoriesWork[r][c];

                            Arena.arrayAttractionWeight[r][c]=stepLengthKernel.vals[looktableR][looktableC]*
                                    (Arena.arrayResourceSelection[r][c]*(familiarity+1));

                    			sum_weights=sum_weights+Arena.arrayAttractionWeight[r][c];
                    		}
                    }

                    if(sum_weights<=0) // only use movement kernel
                    {
                        for(int r=minR;r<maxR;r++)
                        {
                        		looktableR=r-minR;
                        		for(int c=minC;c<maxC;c++)
                        		{
                        			looktableC=c-minC;
                        			Arena.arrayAttractionWeight[r][c]=stepLengthKernel.vals[looktableR][looktableC];
                        			sum_weights=sum_weights+Arena.arrayAttractionWeight[r][c];
                        		}
                        }
                    }


                    // 3. Calculate random step
                    randNumber=drand48()*sum_weights;
                    sumAttractionVectors=0;

                    for(int r=minR;r<maxR;r++)
            			{
            				for(int c=minC;c<maxC;c++)
            				{
            					if(Arena.arrayAttractionWeight[r][c]>0)
            					{
            						if(randNumber>sumAttractionVectors)
            						{
            							if(randNumber<=(sumAttractionVectors+Arena.arrayAttractionWeight[r][c]))
            							{
            								nextCol=c;
            								nextRow=r;
            							}
            						}
            						sumAttractionVectors=sumAttractionVectors+Arena.arrayAttractionWeight[r][c];
            					}
            				}
            			}


                    // 4. Write simulated point
                    totalCount=totalCount+1;
                    totalTraj.animalId.push_back(TrajMetrics.animalId[ind]);
                    totalTraj.run.push_back(run);
                    totalTraj.col.push_back(nextCol);
                    totalTraj.row.push_back(nextRow);

                }
            }
        }

        // 5. Write file
        std::ofstream simulationFile;
        std::string filePath_output2 =  configArg.outputDirectoryPath + "simulations.csv";
        simulationFile.open(filePath_output2.c_str());
        simulationFile<<"animal_id,run,r_patch,c_patch"<<endl;
        	for(int s=0;s<totalTraj.row.size();s++)
        	{
        		simulationFile<<totalTraj.animalId[s]<<","<<totalTraj.run[s]<<","<<totalTraj.row[s]<<","<<totalTraj.col[s]<< std::endl;
        	}

        simulationFile.close();
    }



    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////


    time_t endTime = time(0);

    if(configArg.writingOututs==true)
    {
            MainLogFile<<"\n\n\n" << "Model fitting finished, run time: "<< endTime-iniTime << " seconds" << endl;
            MainLogFile.close();
    }

    cout << "MODEL FITTING ENDS, run time: "<< endTime-iniTime << " seconds" << endl;

}
