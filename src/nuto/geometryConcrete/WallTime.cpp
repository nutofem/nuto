/*
 * WallTime.cpp
 *
 *  Created on: 17 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/WallTime.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include <ctime>
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
