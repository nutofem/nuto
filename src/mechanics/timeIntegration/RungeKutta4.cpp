#ifdef _OPENMP
#include <omp.h>
#endif

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/RungeKutta4.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKutta4::RungeKutta4(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKutta4::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKutta4::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.8 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKutta4::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 4);
    double s;
    switch (rStage)
    {
    case 0:
        s = 0.;
        break;
    case 1:
        s = 0.5;
        break;
    case 2:
        s = 0.5;
        break;
    case 3:
        s = 1.0;
        break;
    default:
        throw Exception("[NuTo::RungeKutta4::GetStageTimeFactor] rStage>3 not implemented.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKutta4::HasTimeChanged(int rStage) const
{
    assert(rStage < 4);
    bool s;
    switch (rStage)
    {
    case 0:
        s = false; // same as last step from the last iteration
        break;
    case 1:
        s = true;
        break;
    case 2:
        s = false;
        break;
    case 3:
        s = true;
        break;
    default:
        throw Exception("[NuTo::RungeKutta4::HasTimeChanged] rStage>3 not implemented.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKutta4::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 4);
    assert(rWeight.size() == 3);
    switch (rStage)
    {
    case 0:
        break;
    case 1:
        rWeight[0] = 0.5;
        break;
    case 2:
        rWeight[0] = 0.0;
        rWeight[1] = 0.5;
        break;
    case 3:
        rWeight[0] = 0.0;
        rWeight[1] = 0.0;
        rWeight[2] = 1.0;
        break;
    default:
        throw Exception("[NuTo::RungeKutta4::GetStageDerivativeFactor] rStage>3 not implemented.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKutta4::GetStageWeights(int rStage) const
{
    assert(rStage < 4);
    double s;
    switch (rStage)
    {
    case 0:
        s = 1. / 6.;
        break;
    case 1:
        s = 1. / 3.;
        break;
    case 2:
        s = 1. / 3.;
        break;
    case 3:
        s = 1. / 6.;
        break;
    default:
        throw Exception("[NuTo::RungeKutta4::GetStageWeights] rStage>3 not implemented.");
    }
    return s;
}
