#pragma once
#include <string>
#include "base/Exception.h"
#include "base/UniqueId.h"

namespace NuTo
{
class DofType : public NuTo::UniqueId<DofType>
{
public:
    DofType(std::string name, int num)
        : mName(name)
        , mNum(num)
    {
        if (mName.empty())
            throw Exception(__PRETTY_FUNCTION__, "Provide a name!");
        if (mNum <= 0)
            throw Exception(__PRETTY_FUNCTION__, "Number of dofs must be greater than zero.");
    }

    const std::string& GetName() const
    {
        return mName;
    }

    int GetNum() const
    {
        return mNum;
    }

private:
    std::string mName;
    int mNum;
};
} /* NuTo */
