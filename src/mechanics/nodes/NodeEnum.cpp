#include "NodeEnum.h"
#include <map>
#include <boost/algorithm/string.hpp>
#include "base/Exception.h"

namespace NuTo
{
namespace Node
{

std::set<eDof> GetDofSet()
{
    const std::set<eDof> set = {eDof::COORDINATES,       eDof::TEMPERATURE,         eDof::DISPLACEMENTS,
                                eDof::NONLOCALEQSTRAIN,  eDof::WATERVOLUMEFRACTION, eDof::RELATIVEHUMIDITY,
                                eDof::ELECTRICPOTENTIAL, eDof::CRACKPHASEFIELD};
    return set;
}


std::map<eDof, std::string> GetDofMap()
{
    std::map<eDof, std::string> attributeMap = {{eDof::COORDINATES, "COORDINATES"},
                                                {eDof::TEMPERATURE, "TEMPERATURE"},
                                                {eDof::DISPLACEMENTS, "DISPLACEMENTS"},
                                                {eDof::NONLOCALEQSTRAIN, "NONLOCALEQSTRAIN"},
                                                {eDof::WATERVOLUMEFRACTION, "WATERVOLUMEFRACTION"},
                                                {eDof::RELATIVEHUMIDITY, "RELATIVEHUMIDITY"},
                                                {eDof::ELECTRICPOTENTIAL, "ELECTRICPOTENTIAL"},
                                                {eDof::CRACKPHASEFIELD, "CRACKPHASEFIELD"}};
    return attributeMap;
}


std::string DofToString(eDof rDof)
{
    try
    {
        return GetDofMap().at(rDof);
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Enum undefined or not implemented.");
    }
}


eDof DofToEnum(std::string dof)
{
    std::string uppercase = boost::to_upper_copy(dof);

    for (const auto& entry : GetDofMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::Exception(__PRETTY_FUNCTION__,
                                   "DofType " + dof + " has no enum equivalent or is not implemented.");
}


bool IsScalar(eDof dofType)
{
    switch (dofType)
    {
    case eDof::COORDINATES:
    case eDof::DISPLACEMENTS:
        return false;
    case eDof::TEMPERATURE:
    case eDof::NONLOCALEQSTRAIN:
    case eDof::RELATIVEHUMIDITY:
    case eDof::WATERVOLUMEFRACTION:
    case eDof::CRACKPHASEFIELD:
    case eDof::ELECTRICPOTENTIAL:
        return true;
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "DOF type not found.");
    }
}


int GetNumComponents(eDof dofType, int dimension)
{
    return IsScalar(dofType) ? 1 : dimension;
}

} // namespace Node
} // namespace NuTo
