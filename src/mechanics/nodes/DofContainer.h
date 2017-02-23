#pragma once

#include <vector>
#include "mechanics/nodes/DofType.h"

namespace NuTo
{
    

template <typename T>
class DofContainer
{
public:
    DofContainer()
    {
        constexpr int magicVectorSize = 10;
        mData.resize(magicVectorSize);
    }
    
    T& operator[](const DofType& rDofType)
    {
        return mData[rDofType.GetId()];
    }
    
    const T& operator[](const DofType& rDofType) const
    {
        return mData[rDofType.GetId()];
    }

private:
    std::vector<T> mData;
};

} /* NuTo */
