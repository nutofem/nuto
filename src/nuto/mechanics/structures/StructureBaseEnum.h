// $Id$
#pragma once
#include <map>
#include <algorithm>

namespace NuTo
{
namespace StructureEnum
{

enum eOutput
{
    HESSIAN0,
    HESSIAN1,
    HESSIAN2,
    HESSIAN2_LUMPED,
    INTERNAL_GRADIENT,
    UPDATE_STATIC_DATA
};

static inline std::map<eOutput, std::string> GetOutputMap()
{
    std::map<eOutput, std::string> map;
    map[HESSIAN0]           = "HESSIAN0";
    map[HESSIAN1]           = "HESSIAN1";
    map[HESSIAN2]           = "HESSIAN2";
    map[HESSIAN2_LUMPED]    = "HESSIAN2_LUMPED";
    map[INTERNAL_GRADIENT]  = "INTERNAL_GRADIENT";
    map[UPDATE_STATIC_DATA] = "UPDATE_STATIC_DATA";
    return map;
}

static inline std::string OutputToString(eOutput rOutput)
{
    try
    {
        return GetOutputMap().find(rOutput)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
    }
}

static inline eOutput OutputToEnum(std::string rOutput)
{
    std::transform(rOutput.begin(), rOutput.end(),rOutput.begin(), ::toupper);

    for(auto entry : GetOutputMap())
        if (entry.second == rOutput)
            return entry.first;

    throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
}

}
}

