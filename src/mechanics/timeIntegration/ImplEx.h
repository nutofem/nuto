/*
 * ImplEx.h
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/timeIntegration/ImplicitExplicitBase.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

namespace NuTo
{

class ImplExCallback;

class ImplEx: public ImplicitExplicitBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#ifndef SWIG
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
    ImplEx() = default;
#endif// SWIG

#endif  // ENABLE_SERIALIZATION
public:

    ImplEx(StructureBase* rStructure);


    virtual ~ImplEx();

    double GetExtrapolationErrorThreshold() const
    {
        return mExtrapolationErrorThreshold;
    }

    void SetExtrapolationErrorThreshold(double rExtrapolationErrorThreshold)
    {
        mExtrapolationErrorThreshold = rExtrapolationErrorThreshold;
    }

    void SetImplExCallback(ImplExCallback* r);

protected:

    //! @brief ... assess the solution and return the new time step
    //! @return ... bool : true - accept solution, false - reject solution
    bool CheckExtrapolationAndAdjustTimeStep() override;

    //! @brief ... extrapolates the static data
    //! @param rTimeStep ... time step object
    void ExtrapolateStaticData(const ConstitutiveTimeStep& rTimeStep) override;

private:

    double mExtrapolationErrorThreshold;

    ImplExCallback* mImplExCallback;

};

} /* namespace NuTo */

