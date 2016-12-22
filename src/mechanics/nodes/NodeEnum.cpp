#include "NodeEnum.h"
#include <boost/algorithm/string.hpp>
#include "mechanics/MechanicsException.h"

const std::set<NuTo::Node::eDof> NuTo::Node::GetDofSet()
{
    const std::set<eDof> set =
    {eDof::COORDINATES,
     eDof::ROTATIONS,
     eDof::TEMPERATURE,
     eDof::DISPLACEMENTS,
     eDof::FINESCALEDISPLACEMENTS,
     eDof::NONLOCALDATA,
     eDof::NONLOCALEQPLASTICSTRAIN,
     eDof::NONLOCALTOTALSTRAIN,
     eDof::NONLOCALEQSTRAIN,
     eDof::WATERVOLUMEFRACTION,
     eDof::RELATIVEHUMIDITY,
     eDof::CRACKPHASEFIELD};
    return set;
}



const std::map<NuTo::Node::eDof, std::string> NuTo::Node::GetDofMap()
{
    const std::map<eDof, std::string> attributeMap =
       {{eDof::COORDINATES,             "COORDINATES"},
        {eDof::ROTATIONS,               "ROTATIONS"},
        {eDof::TEMPERATURE,             "TEMPERATURE"},
        {eDof::DISPLACEMENTS,           "DISPLACEMENTS"},
        {eDof::FINESCALEDISPLACEMENTS,  "FINESCALEDISPLACEMENTS"},
        {eDof::NONLOCALDATA,            "NONLOCALDATA"},
        {eDof::NONLOCALEQPLASTICSTRAIN, "NONLOCALEQPLASTICSTRAIN"},
        {eDof::NONLOCALTOTALSTRAIN,     "NONLOCALTOTALSTRAIN"},
        {eDof::NONLOCALEQSTRAIN,        "NONLOCALEQSTRAIN"},
        {eDof::WATERVOLUMEFRACTION,     "WATERVOLUMEFRACTION"},
        {eDof::RELATIVEHUMIDITY,        "RELATIVEHUMIDITY"},
        {eDof::CRACKPHASEFIELD,         "CRACKPHASEFIELD"}};
    return attributeMap;
}

const std::string NuTo::Node::DofToString(NuTo::Node::eDof rDof)
{
    try
    {
        return GetDofMap().find(rDof)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException("[NuTo::Node::DofToString] Enum undefined or not implemented.");
    }
}

NuTo::Node::eDof NuTo::Node::DofToEnum(std::string rDof)
{
    std::string uppercase = boost::to_upper_copy(rDof);

    for(auto entry : GetDofMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Node::DofToEnum] DofType " + rDof + " has no enum equivalent or is not implemented.");
}
