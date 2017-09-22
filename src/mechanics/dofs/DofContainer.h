#pragma once

#include <map>
#include "mechanics/dofs/DofType.h"

namespace NuTo
{
template <typename T>
class DofContainer
{
public:
    virtual ~DofContainer() = default;

    DofContainer() = default;
    DofContainer(const DofContainer&) = default;
    DofContainer(DofContainer&&) = default;
    DofContainer& operator=(const DofContainer&) = default;
    DofContainer& operator=(DofContainer&&) = default;

    T& operator[](const DofType& dofType)
    {
        return mData[dofType.Id()];
    }

    const T& operator[](const DofType& dofType) const
    {
        return mData.at(dofType.Id());
    }

protected:
    std::map<int, T> mData;
};
} /* NuTo */
