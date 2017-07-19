#ifdef _OPENMP
#include <omp.h>
#endif

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/NystroemQinZhu.h"


//! @brief constructor
//! @param mDimension number of nodes
NuTo::NystroemQinZhu::NystroemQinZhu(StructureBase* rStructure)
    : NystroemBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NystroemQinZhu::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::NystroemQinZhu::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 3);
    double s;
    switch (rStage)
    {
    case 0:
        s = (3. + sqrt(3.)) / 6.;
        break;
    case 1:
        s = (3. - sqrt(3.)) / 6.;
        break;
    case 2:
        s = (3. + sqrt(3.)) / 6.;
        break;
    default:
        throw Exception("[NuTo::NystroemQinZhu::GetStageTimeFactor] error with stage number.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::NystroemQinZhu::HasTimeChanged(int rStage) const
{
    assert(rStage < 3);
    bool s;
    switch (rStage)
    {
    case 0:
        s = true;
        break;
    case 1:
        s = true;
        break;
    case 2:
        s = true;
        break;
    default:
        throw Exception("[NuTo::NystroemQinZhu::HasTimeChanged] error with stage number.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::NystroemQinZhu::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 3);
    assert(rWeight.size() == 2);
    switch (rStage)
    {
    case 0:
        break;
    case 1:
        rWeight[0] = (2. - sqrt(3.)) / 12.;
        break;
    case 2:
        rWeight[0] = 0.0;
        rWeight[1] = sqrt(3.) / 6.;
        break;
    default:
        throw Exception("[NuTo::NystroemQinZhu::GetStageDerivativeFactor] error with stage number.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::NystroemQinZhu::GetStageWeights1(int rStage) const
{
    assert(rStage < 3);
    double s;
    switch (rStage)
    {
    case 0:
        s = (5. - 3. * sqrt(3.)) / 24.;
        break;
    case 1:
        s = (3. + sqrt(3.)) / 12.;
        break;
    case 2:
        s = (1. + sqrt(3.)) / 24.;
        break;
    default:
        throw Exception("[NuTo::NystroemQinZhu::GetStageWeights1] error with stage number.");
    }
    return s;
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::NystroemQinZhu::GetStageWeights2(int rStage) const
{
    assert(rStage < 3);
    double s;
    switch (rStage)
    {
    case 0:
        s = (3. - 2. * sqrt(3.)) / 12.;
        break;
    case 1:
        s = 0.5;
        break;
    case 2:
        s = (3. + 2. * sqrt(3.)) / 12.;
        break;
    default:
        throw Exception("[NuTo::NystroemQinZhu::GetStageWeights2] error with stage number.");
    }
    return s;
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
//! check paper "Lax-Wendroff and Nystrom methods for seismic modelling" by Jing-Bo Chen
double NuTo::NystroemQinZhu::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.58651889452 / std::sqrt(maxGlobalEigenValue);
}
