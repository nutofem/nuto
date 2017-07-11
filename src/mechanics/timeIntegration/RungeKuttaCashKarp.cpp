// $Id: RungeKuttaCashKarp.cpp 575 2011-09-20 18:05:35Z unger3 $

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/RungeKuttaCashKarp.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"


#define orderCashKarp 5

//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKuttaCashKarp::RungeKuttaCashKarp(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKuttaCashKarp::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines (this is wrong do)
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKuttaCashKarp::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.8 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKuttaCashKarp::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 6);
    double s(0);
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
        s = 0.6;
        break;
    case 4:
        s = 1.0;
        break;
    case 5:
        s = 7. / 8.;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKuttaCashKarp::GetStageTimeFactor] rStage<6.");
    }
    return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKuttaCashKarp::HasTimeChanged(int rStage) const
{
    assert(rStage < 6);
    bool s(0);
    switch (rStage)
    {
    case 0:
        s = true; // same as last step from the last iteration
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
    default:
        throw MechanicsException("[NuTo::RungeKuttaCashKarp::HasTimeChanged] rStage<6.");
    }
    return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKuttaCashKarp::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 6);
    assert(rWeight.size() == 5);
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
        rWeight[0] = 0.3;
        rWeight[1] = -0.9;
        rWeight[2] = 1.2;
        break;
    case 4:
        rWeight[0] = -11. / 54.;
        rWeight[1] = 2.5;
        rWeight[2] = -70. / 27.;
        rWeight[3] = 35. / 27.;
        break;
    case 5:
        rWeight[0] = 1631. / 55296.;
        rWeight[1] = 175. / 512.;
        rWeight[2] = 575. / 13824.;
        rWeight[3] = 44275. / 110592.;
        rWeight[4] = 253. / 4096.;
        break;
    default:
        throw MechanicsException("[NuTo::RungeKuttaCashKarp::GetStageDerivativeFactor] rStage<6.");
    }
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKuttaCashKarp::GetStageWeights(int rStage) const
{
    assert(rStage < 6);
    double s;
    if (orderCashKarp == 1)
    {
        switch (rStage)
        {
        case 0:
            s = 1.;
            break;
        case 1:
            s = 0.;
            break;
        case 2:
            s = 0.;
            break;
        case 3:
            s = 0.;
            break;
        case 4:
            s = 0.;
            break;
        case 5:
            s = 0.;
            break;
        default:
            throw MechanicsException("[NuTo::RungeKuttaCashKarp::GetStageWeights] rStage<6.");
        }
    }
    if (orderCashKarp == 2)
    {
        switch (rStage)
        {
        case 0:
            s = -1.5;
            break;
        case 1:
            s = 2.5;
            break;
        case 2:
            s = 0.;
            break;
        case 3:
            s = 0.;
            break;
        case 4:
            s = 0.;
            break;
        case 5:
            s = 0.;
            break;
        default:
            throw MechanicsException("[NuTo::RungeKuttaCashKarp::GetStageWeights] rStage<6.");
        }
    }
    if (orderCashKarp == 3)
    {
        switch (rStage)
        {
        case 0:
            s = 19. / 54.;
            break;
        case 1:
            s = 0.;
            break;
        case 2:
            s = -10. / 27.;
            break;
        case 3:
            s = 55. / 54.;
            break;
        case 4:
            s = 0.;
            break;
        case 5:
            s = 0.;
            break;
        default:
            throw MechanicsException("[NuTo::RungeKuttaCashKarp::GetStageWeights] rStage<6.");
        }
    }
    if (orderCashKarp == 4)
    {
        switch (rStage)
        {
        case 0:
            s = 2825. / 27648.;
            break;
        case 1:
            s = 0.;
            break;
        case 2:
            s = 18575. / 48384.;
            break;
        case 3:
            s = 13525. / 55296.;
            break;
        case 4:
            s = 277. / 14336.;
            break;
        case 5:
            s = 0.25;
            break;
        default:
            throw MechanicsException("[NuTo::RungeKuttaCashKarp::GetStageWeights] rStage<6.");
        }
    }
    if (orderCashKarp == 5)
    {
        switch (rStage)
        {
        case 0:
            s = 37. / 378.;
            break;
        case 1:
            s = 0.;
            break;
        case 2:
            s = 250. / 621.;
            break;
        case 3:
            s = 125. / 594.;
            break;
        case 4:
            s = 0.;
            break;
        case 5:
            s = 512. / 1771.;
            break;
        default:
            throw MechanicsException("[NuTo::RungeKuttaCashKarp::GetStageWeights] rStage<6.");
        }
    }
    return s;
}
