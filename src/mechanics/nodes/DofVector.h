#pragma once

#include "mechanics/nodes/DofContainer.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{

//! @brief dof container that is also capable of performing calculations.
//! @remark We should certainly (!) discuss the name. But I wanted it to be short and catchy :)
class DofVector : public DofContainer<Eigen::VectorXd>
{
public:
    DofVector& operator+=(const DofVector& rOther)
    {
        for (auto i = 0; i < rOther.mData.size(); ++i)
        {
            if (this->mData[i].rows() == 0)
                this->mData[i] = rOther.mData[i];
            else
                this->mData[i] += rOther.mData[i];
        }
        return *this;
    }

    DofVector& operator*=(double rScalar)
    {
        for (auto i = 0; i < this->mData.size(); ++i)
            this->mData[i] *= rScalar;
        return *this;
    }


    friend DofVector operator+(DofVector rLhs, const DofVector& rRhs)
    {
        return std::move(rLhs += rRhs);
    }

    friend DofVector operator*(DofVector rLhs, double rScalar)
    {
        return std::move(rLhs *= rScalar);
    }
};
} /* NuTo */
