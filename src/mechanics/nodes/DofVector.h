#pragma once

#include "mechanics/nodes/DofContainer.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{

//! @brief dof container that is also capable of performing calculations.
template <typename T>
class DofVector : public DofContainer<Eigen::Matrix<T, Eigen::Dynamic, 1>>
{
public:
    DofVector& operator+=(const DofVector& rhs)
    {
        for (const auto& it : rhs.mData)
        {
            auto& thisData = this->mData[it.first];
            if (thisData.rows() == 0)
            {
                thisData = it.second;
            }
            else
            {
                assert(thisData.rows() == it.second.rows());
                thisData += it.second;
            }
        }
        return *this;
    }

    DofVector& operator*=(double scalar)
    {
        for (auto& it : this->mData)
            it.second *= scalar;
        return *this;
    }

    friend DofVector operator+(DofVector lhs, const DofVector& rhs)
    {
        return std::move(lhs += rhs);
    }

    friend DofVector operator*(DofVector lhs, double scalar)
    {
        return std::move(lhs *= scalar);
    }
    
    friend std::ostream& operator<<(std::ostream& out, const DofVector<T>& dofVector)
    {
        for (auto const& data : dofVector.mData)
        {
            out << "====" << std::endl;
            out << data.second << std::endl;
        }
        out << "====" << std::endl;
        return out;
    }
};
} /* NuTo */
