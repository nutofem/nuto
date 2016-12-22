#include "StructureBaseEnum.h"

#include <algorithm>
#include "mechanics/MechanicsException.h"

const std::map<NuTo::eStructureOutput, std::string> NuTo::GetOutputMap()
{
    const std::map<eStructureOutput, std::string> map =
       {{eStructureOutput::HESSIAN0,"HESSIAN0"},
        {eStructureOutput::HESSIAN1,"HESSIAN1"},
        {eStructureOutput::HESSIAN2,"HESSIAN2"},
        {eStructureOutput::HESSIAN2_LUMPED,"HESSIAN2_LUMPED"},
        {eStructureOutput::INTERNAL_GRADIENT,"INTERNAL_GRADIENT"},
        {eStructureOutput::UPDATE_STATIC_DATA, "UPDATE_STATIC_DATA"}};
    return map;
}



std::string NuTo::StructureOutputToString(NuTo::eStructureOutput rOutput)
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



NuTo::eStructureOutput NuTo::StructureOutputToEnum(std::string rOutput)
{
    std::transform(rOutput.begin(), rOutput.end(),rOutput.begin(), ::toupper);

    for(auto entry : GetOutputMap())
        if (entry.second == rOutput)
            return entry.first;

    throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
}
