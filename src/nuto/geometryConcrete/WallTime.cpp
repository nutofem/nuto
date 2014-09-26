/*
 * WallTime.cpp
 *
 *  Created on: 17 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/WallTime.h"
#include <ctime>


#ifdef _OPENMP
#include <omp.h>
#endif

//! @return ... current wall time, parallel or sequential
double NuTo::WallTime::Get()
{
#ifdef _OPENMP
	return omp_get_wtime();
#else
	return static_cast<double>(clock()) / CLOCKS_PER_SEC;
#endif
}

std::string NuTo::WallTime::TimeStamp()
{
    char buffer[120];
    time_t now = time(0);
    tm * timeinfo = localtime(&now);
    strftime(buffer, 120, "%F_%H-%M-%S", timeinfo);
    return std::string(buffer);
}
