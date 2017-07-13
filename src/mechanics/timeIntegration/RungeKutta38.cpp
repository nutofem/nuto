// $Id: RungeKutta38.cpp 575 2011-09-20 18:05:35Z unger3 $

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/RungeKutta38.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"


//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKutta38::RungeKutta38(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKutta38::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKutta38::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.8 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKutta38::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 4);
    double s(0);
    switch (rStage)
    {
    case 0:
        s = 0.;
        break;
    case 1:
        s = 1. / 3.;
        break;
    case 2:
        s = 2. / 3.;
        break;
    case 3:
        s = 1.0;
        break;
    default:
        throw Exception("[NuTo::RungeKutta38::GetStageTimeFactor] rStage<4.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKutta38::HasTimeChanged(int rStage) const
{
    assert(rStage < 4);
    bool s(0);
    switch (rStage)
    {
    case 0:
        s = false; // same as last step from the last iteration
        break;
    case 1:
        s = true;
        break;
    case 2:
        s = true;
        break;
    case 3:
        s = true;
        break;
    default:
        throw Exception("[NuTo::RungeKutta38::HasTimeChanged] rStage<4.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKutta38::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 4);
    assert(rWeight.size() == 3);
    switch (rStage)
    {
    case 0:
        break;
    case 1:
        rWeight[0] = 1. / 3.;
        break;
    case 2:
        rWeight[0] = -1. / 3.;
        rWeight[1] = 1.;
        break;
    case 3:
        rWeight[0] = 1.0;
        rWeight[1] = -1.0;
        rWeight[2] = 1.0;
        break;
    default:
        throw Exception("[NuTo::RungeKutta38::GetStageDerivativeFactor] rStage<4.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKutta38::GetStageWeights(int rStage) const
{
    assert(rStage < 4);
    double s(0);
    switch (rStage)
    {
    case 0:
        s = 1. / 8.;
        break;
    case 1:
        s = 3. / 8.;
        break;
    case 2:
        s = 3. / 8.;
        break;
    case 3:
        s = 1. / 8.;
        break;
    default:
        throw Exception("[NuTo::RungeKutta38::GetStageWeights] rStage<4.");
    }
    return s;
}
