#pragma once

#include <map>
#include <eigen3/Eigen/Core>
#include "mechanics/dofs/DofType.h"

namespace NuTo
{

template <typename T>
class DofVector
{
public:
    //! @brief const access
    //! @param dofType dof type
    //! @return const reference to existing value, throws if there is no value
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& operator[](const DofType& dofType) const
    {
        return mData.at(dofType);
    }

    //! @brief nonconst access
    //! @param dofType dof type
    //! @return reference to either
    //!         an existing value
    //!         an newly default constructed value
    //! @remark This requires T to be default constructable.
    Eigen::Matrix<T, Eigen::Dynamic, 1>& operator[](const DofType& dofType)
    {
        return mData[dofType];
    }

    //! @brief performs _uninitialized addition_ that resizes the data to the length of `rhs`
    DofVector& operator+=(const DofVector& rhs)
    {
        for (auto& entry : rhs.mData)
        {
            if (mData.find(entry.first) != mData.end())
                mData[entry.first] += entry.second;
            else
                mData[entry.first] = entry.second;
        }
        return *this;
    }

    DofVector& operator*=(double scalar)
    {
        for (auto& entry : mData)
            entry.second *= scalar;
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

    friend std::ostream& operator<<(std::ostream& out, const DofVector& v)
    {
        for (auto& entry : v.mData)
        {
            out << "==== " << entry.first.GetName() << "====\n";
            out << entry.second << '\n';
        }
        return out;
    }

    //! @brief calculates (*this) += `rhs` * `scalar`
    //! @remark This is the common case in the numerical integration - summing up all integration point contributions
    //! scaled with DetJ and the integration point weight.
    void AddScaled(const DofVector& rhs, double scalar)
    {
        for (auto& entry : rhs.mData)
        {
            if (mData.find(entry.first) != mData.end())
                mData[entry.first] += entry.second * scalar;
            else
                mData[entry.first] = entry.second * scalar;
        }
    }


protected:
    //! @brief data container
    std::map<DofType, Eigen::Matrix<T, Eigen::Dynamic, 1>, CompareDofType>
            mData;
};

} /* NuTo */
