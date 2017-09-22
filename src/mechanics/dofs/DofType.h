#pragma once
#include <string>
#include "base/Exception.h"

namespace NuTo
{
class DofType
{
public:
    DofType(std::string name, int num, int id)
        : mName(name)
        , mNum(num)
        , mId(id)
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

    int GetId() const
    {
        return mId;
    }

private:
    std::string mName;
    int mNum;
    int mId;
};
} /* NuTo */
