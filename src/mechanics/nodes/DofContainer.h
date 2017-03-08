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

    virtual ~DofContainer()           = default;
    DofContainer(const DofContainer&) = default;
    DofContainer(DofContainer&&)      = default;
    DofContainer& operator=(const DofContainer&) = default;
    DofContainer& operator=(DofContainer&&) = default;

    T& operator[](const DofType& rDofType)
    {
        //return mData[rDofType.GetId()];
        return mData.at(rDofType.GetId());
    }

    const T& operator[](const DofType& rDofType) const
    {
        //return mData[rDofType.GetId()];
        return mData.at(rDofType.GetId());
    }

protected:
    std::vector<T> mData;
};

} /* NuTo */
