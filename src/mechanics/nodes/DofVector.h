#pragma once

#include "mechanics/nodes/DofContainer.h"
#include <Eigen/Core>

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
                thisData = it.second;
            else
                thisData += it.second;
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
        lhs += rhs;
        return lhs;
    }

    friend DofVector operator*(DofVector lhs, double scalar)
    {
        lhs *= scalar;
        return lhs;
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
