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

    T& operator[](const DofType& dofType)
    {
        return mData[dofType.Id()];
    }

    const T& operator[](const DofType& dofType) const
    {
        return mData.at(dofType.Id());
    }

    void Insert(const DofType& dofType, T t)
    {
        mData.emplace(dofType.Id(), t);
    }

protected:
    std::map<int, T> mData;
};
} /* NuTo */
