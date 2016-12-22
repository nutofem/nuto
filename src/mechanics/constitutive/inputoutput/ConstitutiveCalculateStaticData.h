/*
 * ConstitutiveCalculateStaticData.h
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"


namespace NuTo
{


enum class eCalculateStaticData
{
    EULER_FORWARD,
    EULER_BACKWARD,
    USE_PREVIOUS
};



class ConstitutiveCalculateStaticData: public ConstitutiveIOBase
{
public:
    ConstitutiveCalculateStaticData(eCalculateStaticData rCalculateStaticData, int rIndexOfPreviousStaticData = 0) :
        ConstitutiveIOBase() ,
        mCalculateStaticData(rCalculateStaticData),
        mIndexOfPreviousStaticData(rIndexOfPreviousStaticData)
    {}

    virtual std::unique_ptr<ConstitutiveIOBase> clone() override
    {
        return std::make_unique<ConstitutiveCalculateStaticData>(*this);
    }

    //! @brief returns (rXn + rTimeStep[0] / rTimeStep[1] * (rXn - rXn_m1)
    template <typename T>
    static T EulerForward(const T& rXn, const T& rXn_m1, const ConstitutiveIOBase& rTimeStep)
    {
        assert (rTimeStep.GetNumRows() >= 2 && "At least two time steps required. Current and previous one.");

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

    eCalculateStaticData GetCalculateStaticData() const
    {
        return mCalculateStaticData;
    }

    void SetCalculateStaticData(const eCalculateStaticData& rCalculateStaticData)
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
    eCalculateStaticData mCalculateStaticData;
    int mIndexOfPreviousStaticData;
};

} /* namespace NuTo */

