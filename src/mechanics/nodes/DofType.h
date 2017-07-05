#pragma once
#include <string>
#include "mechanics/MechanicsException.h"

namespace NuTo
{
class DofType
{
public:
    DofType(std::string rName, int rNum, int rId)
        : mName(rName)
        , mNum(rNum)
        , mId(rId)
    {
        if (mName.empty())
            throw MechanicsException(__PRETTY_FUNCTION__, "Provide a name!");
        if (mNum <= 0)
            throw MechanicsException(__PRETTY_FUNCTION__, "Number of dofs must be greater than zero.");
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
