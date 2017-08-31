#pragma once

#include <map>
#include <eigen3/Eigen/Core>
#include "mechanics/nodes/DofType.h"

namespace NuTo
{

//! @brief dof container that is also capable of performing calculations.
template <typename T>
class DofMatrixContainer
{
public:
    DofMatrixContainer& operator+=(const DofMatrixContainer& rhs)
    {
        for (const auto& it : rhs.mData)
        {
            auto& thisData = mData[it.first];
            if (thisData.size() == 0)
                thisData = it.second;
            else
                thisData += it.second;
        }
        return *this;
    }

    DofMatrixContainer& operator*=(double scalar)
    {
        for (auto& it : mData)
            it.second *= scalar;
        return *this;
    }

    T& operator()(const DofType& d0, const DofType& d1)
    {
        return mData[CantorParingFunction(d0.GetId(), d1.GetId())];
    }

    const T& operator()(const DofType& d0, const DofType& d1) const
    {
        return mData.at(CantorParingFunction(d0.GetId(), d1.GetId()));
    }

    friend DofMatrixContainer operator+(DofMatrixContainer lhs, const DofMatrixContainer& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    friend DofMatrixContainer operator*(DofMatrixContainer lhs, double scalar)
    {
        lhs *= scalar;
        return lhs;
    }

    friend std::ostream& operator<<(std::ostream& out, const DofMatrixContainer<T>& dofMatrix)
    {
        for (auto const& data : dofMatrix.mData)
        {
            auto xy = CantorPairingFunctionReverse(data.first);
            out << "=== " << xy.first << " " << xy.second << " ==="  << std::endl;
            out << data.second << std::endl;
        }
        out << "====" << std::endl;
        return out;
    }

private:
    //! @brief https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
    //! maps NxN --> N uniquely, deterministic, nice.
    static int CantorParingFunction(int k1, int k2)
    {
        return k2 + (k1 + k2) * (k1 + k2 + 1) / 2;
    }

    static std::pair<int, int> CantorPairingFunctionReverse(int z)
    {
        int w = std::floor((std::sqrt(8 * z + 1) - 1) / 2);
        int t = (w * w + w) / 2;
        return std::make_pair(z - t, w - z + t);
    }

    std::map<int, T> mData;
};
} /* NuTo */
