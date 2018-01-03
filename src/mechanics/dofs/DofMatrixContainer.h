#pragma once

#include "mechanics/dofs/DofCalcContainer.h"

namespace NuTo
{

//! @brief dof container that is also capable of performing calculations.
template <typename T>
class DofMatrixContainer : public DofCalcContainer<T>
{
public:
    T& operator()(DofType d0, DofType d1)
    {
        return this->ResizingIdAccess(CantorParingFunction(d0.Id(), d1.Id()));
    }

    const T& operator()(DofType d0, DofType d1) const
    {
        return this->mData[CantorParingFunction(d0.Id(), d1.Id())];
    }

    const T& operator[](DofType) const
    {
        throw Exception(__PRETTY_FUNCTION__, "This method is for vectors only, not for matrices.");
    }

    friend std::ostream& operator<<(std::ostream& out, const DofMatrixContainer<T>& dofMatrix)
    {
        for (size_t i = 0; i < dofMatrix.mData.size(); ++i)
        {
            auto xy = CantorPairingFunctionReverse(i);
            out << "=== " << xy.first << " " << xy.second << " ===" << std::endl;
            out << dofMatrix.mData[i] << std::endl;
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
};
} /* NuTo */
