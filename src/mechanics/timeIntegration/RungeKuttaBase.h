
#pragma once

#include "mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (Newmark, but you can use it for statics as well with setting
//! the flag isDynamic to false)
class RungeKuttaBase : public TimeIntegrationBase
{

public:
    //! @brief constructor
    RungeKuttaBase(StructureBase* rStructure);

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    void Solve(double rTimeDelta) override;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override;

    //! @brief ... return number of intermediate stages
    virtual int GetNumStages() const = 0;

    //! @brief ... return delta time factor of intermediate stages (c in Butcher tableau)
    virtual double GetStageTimeFactor(int rStage) const = 0;

    //! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
    virtual void GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const = 0;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    virtual double GetStageWeights(int rStage) const = 0;

    //! @brief ... return, if time (e.g. for the calculation of external loads) has changed
    virtual bool HasTimeChanged(int rStage) const = 0;


    // temporary fix
    //! @brief sets the  time step for the time integration procedure (initial value)
    virtual void SetTimeStep(double rTimeStep) override
    {
        mTimeStep = rTimeStep;
    }

    //! @brief returns the  time step for the time integration procedure (current value)
    virtual double GetTimeStep() const override
    {
        return mTimeStep;
    }

protected:
    double mTime = 0.;
    double mTimeStep = 0.;
};
} // namespace NuTo
