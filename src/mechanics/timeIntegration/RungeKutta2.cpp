# ifdef _OPENMP
#include <omp.h>
# endif

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/RungeKutta2.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"


//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKutta2::RungeKutta2 (StructureBase* rStructure)  : RungeKuttaBase (rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKutta2::Info()const
{
	TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKutta2::CalculateCriticalTimeStep()const
{
	double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.0/std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKutta2::GetStageTimeFactor(int rStage)const
{
    assert(rStage<2);
	double s;
	switch(rStage)
	{
	case 0:
		s = 0.;
		break;
	case 1:
        s = 0.5;
		break;
	default:
        throw MechanicsException ( "[NuTo::RungeKutta2::GetStageTimeFactor] rStage>3 not implemented." );
	}
	return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKutta2::HasTimeChanged(int rStage)const
{
    assert(rStage<2);
	bool s;
	switch(rStage)
	{
	case 0:
		s = false; //same as last step from the last iteration
		break;
    case 1:
		s = true;
		break;
	default:
        throw MechanicsException ( "[NuTo::RungeKutta2::HasTimeChanged] rStage>3 not implemented." );
	}
	return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKutta2::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const
{
    assert(rStage<2);
    assert(rWeight.size()==1);
	switch(rStage)
	{
	case 0:
		break;
	case 1:
        rWeight[0] = 0.5;
		break;
	default:
        throw MechanicsException ( "[NuTo::RungeKutta2::GetStageDerivativeFactor] rStage>3 not implemented." );
	}
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKutta2::GetStageWeights(int rStage)const
{
    assert(rStage<2);
	double s;
	switch(rStage)
	{
	case 0:
        s = 0.0;
		break;
	case 1:
        s = 1.0;
		break;
	default:
        throw MechanicsException ( "[NuTo::RungeKutta2::GetStageWeights] rStage>3 not implemented." );
	}
	return s;
}
