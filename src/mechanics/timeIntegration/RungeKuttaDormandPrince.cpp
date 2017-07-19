#ifdef _OPENMP
#include <omp.h>
#endif

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/RungeKuttaDormandPrince.h"


//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKuttaDormandPrince::RungeKuttaDormandPrince(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKuttaDormandPrince::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines (this is wrong do)
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKuttaDormandPrince::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.8 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKuttaDormandPrince::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 7);
    double s;
    switch (rStage)
    {
    case 0:
        s = 0.;
        break;
    case 1:
        s = 0.2;
        break;
    case 2:
        s = 0.3;
        break;
    case 3:
        s = 0.8;
        break;
    case 4:
        s = 8. / 9.;
        break;
    case 5:
        s = 1.0;
        break;
    case 6:
        s = 1.0;
        break;
    default:
        throw Exception("[NuTo::RungeKuttaDormandPrince::GetStageTimeFactor] rStage<7.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKuttaDormandPrince::HasTimeChanged(int rStage) const
{
    assert(rStage < 7);
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
        s = true;
        break;
    case 3:
        s = true;
        break;
    case 4:
        s = true;
        break;
    case 5:
        s = true;
        break;
    case 6:
        s = false;
        break;
    default:
        throw Exception("[NuTo::RungeKuttaDormandPrince::HasTimeChanged] rStage<7.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKuttaDormandPrince::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 7);
    assert(rWeight.size() == 6);
    switch (rStage)
    {
    case 0:
        break;
    case 1:
        rWeight[0] = 0.2;
        break;
    case 2:
        rWeight[0] = 3. / 40.;
        rWeight[1] = 9. / 40.;
        break;
    case 3:
        rWeight[0] = 44. / 45.;
        rWeight[1] = -56. / 15.;
        rWeight[2] = 32. / 9.;
        break;
    case 4:
        rWeight[0] = 19372. / 6561.;
        rWeight[1] = -25360. / 2187.;
        rWeight[2] = 64448. / 6561.;
        rWeight[3] = -212. / 729.;
        break;
    case 5:
        rWeight[0] = 9017. / 3168.;
        rWeight[1] = -355. / 33.;
        rWeight[2] = 46732. / 5247.;
        rWeight[3] = 49. / 176.;
        rWeight[4] = -5103. / 18656.;
        break;
    case 6:
        rWeight[0] = 35. / 384.;
        rWeight[1] = 0.;
        rWeight[2] = 500. / 1113.;
        rWeight[3] = 125. / 192.;
        rWeight[4] = -2187. / 6784.;
        rWeight[5] = 11. / 84.;
        break;
    default:
        throw Exception("[NuTo::RungeKuttaDormandPrince::GetStageDerivativeFactor] rStage<7.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKuttaDormandPrince::GetStageWeights(int rStage) const
{
    assert(rStage < 7);
    double s;
    switch (rStage)
    {
    case 0:
        s = 35. / 385;
        break;
    case 1:
        s = 0.;
        break;
    case 2:
        s = 500. / 1113.;
        break;
    case 3:
        s = 125. / 192.;
        break;
    case 4:
        s = -2187. / 6784.;
        break;
    case 5:
        s = 11. / 84.;
        break;
    case 6:
        s = 0.;
        break;
    default:
        throw Exception("[NuTo::RungeKuttaDormandPrince::GetStageWeights] rStage<7.");
    }
    return s;
}
