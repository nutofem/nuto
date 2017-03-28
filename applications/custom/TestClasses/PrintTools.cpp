#include "PrintTools.h"

#include <iostream>


void PrintTools::printArray_int(int *rArray, int rSize, std::__cxx11::string rTitle, int rPID)
{
    std::cout << rTitle << "\n" << std::endl;
    for (int j = 0; j < rSize; ++j)
    {
        std::cout << "PID: " << rPID << ":  " << j << " --> " << rArray[j]  << std::endl;
    }
}


void PrintTools::printVector(std::vector<int> rVector, std::__cxx11::string rTitle, int rPID)
{
    std::cout << rTitle << "\n" << std::endl;
    for (int j = 0; j < rVector.size(); ++j)
    {
        std::cout << "PID: " << rPID << ":  " << j << " --> " << rVector[j]  << std::endl;
    }
}


void PrintTools::printMap_int_int(std::map<int, int> rMap, std::__cxx11::string rTitle, int rPID)
{
    std::cout << rTitle << "\n" << std::endl;
    for (int j = 0; j < rMap.size(); ++j)
    {
        std::cout << "PID: " << rPID << ":  " << j << " --> [" << rMap[j] << "]" << std::endl;
    }
}

