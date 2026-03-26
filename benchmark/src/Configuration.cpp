#include "header.h"


// Reads the configuration file
structConfig launchConfig_txt_sep(std::string path_to_config)
{
    structConfig returnValues;
    std::string directory_rasters, filePath, line, extracted;
    std::string filePath_sLength = path_to_config;
    std::string line2, extracted2, line3;
    std::ifstream sLength2 (filePath_sLength.c_str());
    getline(sLength2, line2);


    // GET THE DIRECTORIES
    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.outputDirectoryPath=line2;

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.rasterDirectoryPath=line2;

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.trajectoryPath=line2;


   // GET THE WRITING OPTIONS
    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    if(line2=="true"){returnValues.writingOututs=true;}
    if(line2=="false"){returnValues.writingOututs=false;}


    // GET THE SCOPE
    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.modelType=line2;


    // GET THE APPROXIMATION
    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.thresholdApproxKernel=std::stoi(line2);

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.thresholdMemoryKernel=std::stoi(line2);


    // GET RASTER RESOURCES (variable length) + ASSOCIATED SELECTION COEFFICIENTS
    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);        // remove the first column
    size_t n = std::count(line2.begin(), line2.end(), ',')+1;
    int i;
    for(i=1;i<=n;i++)
    {
        if(i<n)
        {
            extracted=line2.substr(0,line2.find(",",0));
            returnValues.resourceNames.push_back(extracted);
            line2=line2.substr(line2.find(",",0)+1,line2.length()-1);
        }
        else
        {
            extracted=line2;
            returnValues.resourceNames.push_back(line2);
        }
    }

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    for(i=1;i<=n;i++)
    {
        if(i<n)
        {
            extracted=line2.substr(0,line2.find(",",0));
            returnValues.selectionCoef.push_back(std::stod(extracted));
            line2=line2.substr(line2.find(",",0)+1,line2.length()-1);
        }
        else
        {
            extracted=line2;
            returnValues.selectionCoef.push_back(std::stod(extracted));
        }
    }


    // GET MEMORY PARAMETERS
    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.memoryRL=std::stod(line2);

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.memoryWL=std::stod(line2);

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.memoryRD=std::stod(line2);
    returnValues.memoryRD_cplm=(1-returnValues.memoryRD);

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.memoryWD=std::stod(line2);
    returnValues.memoryWD_cplm=(1-returnValues.memoryWD);


    // GET DISTANCE PARAMETERS
    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.memoryRDist=std::stod(line2);
    returnValues.memoryRDist=returnValues.memoryRDist*(-1);

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.memoryWDist=std::stod(line2);
    returnValues.memoryWDist=returnValues.memoryWDist*(-1);

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.stepLengthDist=std::stod(line2);

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.stepLengthShape=std::stod(line2);


    // OTHER PARAMETERS
    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.nSimulatedSteps=std::stoi(line2);

    getline(sLength2, line2);
    line2=line2.erase(0, line2.find(";")+1);
    returnValues.nSimulatedRuns=std::stoi(line2);


    return returnValues;
}
