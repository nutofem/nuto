/*
 * ImplEx.h
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/timeIntegration/ImplicitExplicitBase.h"

namespace NuTo
{

class ImplExCallback;

class ImplEx : public ImplicitExplicitBase
{
public:
    ImplEx(StructureBase* rStructure);


    virtual ~ImplEx() = default;

    double GetExtrapolationErrorThreshold() const
    {
        return mExtrapolationErrorThreshold;
    }

    void SetExtrapolationErrorThreshold(double rExtrapolationErrorThreshold)
    {
        mExtrapolationErrorThreshold = rExtrapolationErrorThreshold;
    }

    void SetImplExCallback(std::shared_ptr<ImplExCallback> r);


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
    //! @brief ... assess the solution and return the new time step
    //! @return ... bool : true - accept solution, false - reject solution
    bool CheckExtrapolationAndAdjustTimeStep() override;

    //! @brief ... extrapolates the static data
    //! @param rTimeStep ... time step object
    void ExtrapolateStaticData(const ConstitutiveTimeStep& rTimeStep) override;

private:
    double mExtrapolationErrorThreshold;

    std::shared_ptr<ImplExCallback> mImplExCallback;
};

} /* namespace NuTo */
