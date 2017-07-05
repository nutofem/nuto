#pragma once

#include "mechanics/timeIntegration/RungeKuttaBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (Newmark, but you can use it for statics as well with setting the flag isDynamic to false)
class RungeKutta4 : public RungeKuttaBase
{

public:

    //! @brief constructor
    RungeKutta4(StructureBase* rStructure);

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep()const override
    {
    	return true;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep()const override;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const override;

    //! @brief ... return number of intermediate stages
    int GetNumStages()const override
    {
    	return 4;
    }

    //! @brief ... return delta time factor of intermediate stages (c in Butcher tableau)
    double GetStageTimeFactor(int rStage)const override;

    //! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
    void GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const override;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    double GetStageWeights(int rStage)const override;

    //! @brief ... return, if time (e.g. for the calculation of external loads) has changed
    bool HasTimeChanged(int rStage)const override;

protected:
    //empty private construct required for serialization
};
} //namespace NuTo



