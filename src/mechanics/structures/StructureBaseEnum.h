#pragma once
#include <map>
#include <string>


namespace NuTo
{


enum class eStructureOutput
{
    HESSIAN0,
    HESSIAN1,
    HESSIAN2,
    HESSIAN2_LUMPED,
    INTERNAL_GRADIENT,
    UPDATE_STATIC_DATA
};

const std::map<eStructureOutput, std::string> GetOutputMap();
std::string StructureOutputToString(eStructureOutput rOutput);
eStructureOutput StructureOutputToEnum(std::string rOutput);


} // namespace NuTo
