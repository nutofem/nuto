/*
 * ConstitutiveCalculateStaticData.h
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"


namespace NuTo
{

namespace CalculateStaticData
{

    enum eCalculateStaticData
    {
        EULER_FORWARD,
        EULER_BACKWARD,
        USE_PREVIOUS
    };

}  // namespace CalculateStaticData


class ConstitutiveCalculateStaticData: public ConstitutiveIOBase
{
public:
    ConstitutiveCalculateStaticData(CalculateStaticData::eCalculateStaticData rCalculateStaticData, int rIndexOfPreviousStaticData = 0) :
        ConstitutiveIOBase() ,
        mCalculateStaticData(rCalculateStaticData),
        mIndexOfPreviousStaticData(rIndexOfPreviousStaticData)
    {}


    //! @brief returns (rXn + rTimeStep[0] / rTimeStep[1] * (rXn - rXn_m1)
    template <typename T>
    static T EulerForward(const T& rXn, const T& rXn_m1, const ConstitutiveIOBase& rTimeStep)
    {
        rTimeStep.AssertIsVector<2>(Constitutive::Output::UPDATE_STATIC_DATA, __PRETTY_FUNCTION__); // enum is just a dummy... I'm sorry.

        // this is a bit wierd since this method should be used for ConstiutiveScalar/Vector/Matrix,
        // all inheriting from Eigen::Matrix<double>
        // especially the return types of Eigen, aka Eigen::CwiseUnaryOp<Eigen::internal::scalar_multiple_op<double> ...
        // cannot be converted to ConstitutiveIOBase.

        T result = rXn;
        result -= rXn_m1;
        result *= rTimeStep[0] / rTimeStep[1];
        result += rXn;
        return result;
    }

    CalculateStaticData::eCalculateStaticData GetCalculateStaticData() const
    {
        return mCalculateStaticData;
    }

    void SetCalculateStaticData(const CalculateStaticData::eCalculateStaticData& rCalculateStaticData)
    {
        mCalculateStaticData = rCalculateStaticData;
    }

    int GetIndexOfPreviousStaticData() const
    {
        return mIndexOfPreviousStaticData;
    }

    void SetIndexOfPreviousStaticData(int rIndexOfPreviousStaticData)
    {
        mIndexOfPreviousStaticData = rIndexOfPreviousStaticData;
    }

private:
    CalculateStaticData::eCalculateStaticData mCalculateStaticData;
    int mIndexOfPreviousStaticData;
};

} /* namespace NuTo */

