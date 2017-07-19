#ifdef _OPENMP
#include <omp.h>
#endif

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/RungeKutta3.h"


//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKutta3::RungeKutta3(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKutta3::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKutta3::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.51 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKutta3::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 3);
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
        s = 1.;
        break;
    default:
        throw Exception("[NuTo::RungeKutta3::GetStageTimeFactor] rStage>3 not implemented.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKutta3::HasTimeChanged(int rStage) const
{
    assert(rStage < 3);
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
    default:
        throw Exception("[NuTo::RungeKutta3::HasTimeChanged] rStage>3 not implemented.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKutta3::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 3);
    assert(rWeight.size() == 2);
    switch (rStage)
    {
    case 0:
        break;
    case 1:
        rWeight[0] = 0.5;
        break;
    case 2:
        rWeight[0] = -1.0;
        rWeight[1] = 2.0;
        break;
    default:
        throw Exception("[NuTo::RungeKutta3::GetStageDerivativeFactor] rStage>3 not implemented.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKutta3::GetStageWeights(int rStage) const
{
    assert(rStage < 3);
    double s;
    switch (rStage)
    {
    case 0:
        s = 1. / 6.;
        break;
    case 1:
        s = 2. / 3.;
        break;
    case 2:
        s = 1. / 6.;
        break;
    default:
        throw Exception("[NuTo::RungeKutta3::GetStageWeights] rStage>3 not implemented.");
    }
    return s;
}
