#pragma once

#include <vector>
#include <eigen3/Eigen/Core>
#include "mechanics/dofs/DofType.h"

namespace NuTo
{

//! @brief dof container that is also capable of performing calculations.
template <typename T>
class DofCalcContainer
{
public:
    const T& operator[](const DofType& dofType) const
    {
        return mData.at(dofType.Id());
    }

    T& operator[](const DofType& dofType)
    {
        return ResizingIdAccess(dofType.Id());
    }

    DofCalcContainer& operator+=(const DofCalcContainer& rhs)
    {
        for (int i = 0; i < rhs.mData.size(); ++i)
        {
            auto& thisData = ResizingIdAccess(i);
            if (thisData.size() == 0)
                thisData = rhs.mData[i];
            else
                thisData += rhs.mData[i];
        }
        return *this;
    }

    void AddScaled(const DofCalcContainer& rhs, double scalar)
    {
        for (int i = 0; i < rhs.mData.size(); ++i)
        {
            auto& thisData = ResizingIdAccess(i);
            if (thisData.size() == 0)
                thisData = rhs.mData[i] * scalar;
            else
                thisData += rhs.mData[i] * scalar;
        }
    }

    DofCalcContainer& operator*=(double scalar)
    {
        for (auto& data : this->mData)
            data *= scalar;
        return *this;
    }

    friend DofCalcContainer operator+(DofCalcContainer lhs, DofCalcContainer&& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    friend DofCalcContainer operator*(DofCalcContainer lhs, double scalar)
    {
        lhs *= scalar;
        return lhs;
    }

    friend std::ostream& operator<<(std::ostream& out, const DofCalcContainer<T>& v)
    {
        for (int i = 0; i < v.mData.size(); ++i)
        {
            out << "==== Index " << i << "====\n";
            out << v.mData[i] << '\n';
        }
        return out;
    }

protected:
    std::vector<T> mData;

    //! @brief access to `i`-th entry of the data. Default constructs data up to `i`, if needed.
    //! @param i data index
    //! @return nonconst reference to `i`-th data entry
    T& ResizingIdAccess(int i)
    {
        if (mData.size() < i + 1)
            mData.resize(i + 1);
        return mData[i];
    }
};
} /* NuTo */
