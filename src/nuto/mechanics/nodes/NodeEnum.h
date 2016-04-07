// $Id$ 
#ifndef NODEENUM_H_
#define NODEENUM_H_

#include <map>
#include <boost/algorithm/string.hpp>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
namespace Node
{
enum eNodeType
{
    NodeCoordinates,
    NodeCoordinatesDof,
    NodeCoordinatesDofNonlocalData,
};

enum eDof : unsigned char
{
    COORDINATES,
    ROTATIONS,
    TEMPERATURES,
    DISPLACEMENTS,
    FINESCALEDISPLACEMENTS,
    NONLOCALDATA,
    NONLOCALEQPLASTICSTRAIN,
    NONLOCALTOTALSTRAIN,
    NONLOCALEQSTRAIN,
    WATERVOLUMEFRACTION,
    RELATIVEHUMIDITY
};

constexpr size_t maxDoFEnumValue = std::numeric_limits<unsigned char>::max();

constexpr unsigned short CombineDofs(eDof rDof1, eDof rDof2)
{
    return rDof1 * maxDoFEnumValue + rDof2;
}

static inline std::map<eDof, std::string> GetDofMap()
{
    std::map<eDof, std::string> attributeMap;
    attributeMap[eDof::COORDINATES]             = "COORDINATES";
    attributeMap[eDof::ROTATIONS]               = "ROTATIONS";
    attributeMap[eDof::TEMPERATURES]            = "TEMPERATURES";
    attributeMap[eDof::DISPLACEMENTS]           = "DISPLACEMENTS";
    attributeMap[eDof::FINESCALEDISPLACEMENTS]  = "FINESCALEDISPLACEMENTS";
    attributeMap[eDof::NONLOCALDATA]            = "NONLOCALDATA";
    attributeMap[eDof::NONLOCALEQPLASTICSTRAIN] = "NONLOCALEQPLASTICSTRAIN";
    attributeMap[eDof::NONLOCALTOTALSTRAIN]     = "NONLOCALTOTALSTRAIN";
    attributeMap[eDof::NONLOCALEQSTRAIN]        = "NONLOCALEQSTRAIN";
    attributeMap[eDof::WATERVOLUMEFRACTION]     = "WATERVOLUMEFRACTION";
    attributeMap[eDof::RELATIVEHUMIDITY]        = "RELATIVEHUMIDITY";
    return attributeMap;
}


static inline std::string DofToString(eDof rDof)
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

static inline eDof DofToEnum(std::string rDof)
{
    std::string uppercase = boost::to_upper_copy(rDof);

    for(auto entry : GetDofMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Node::DofToEnum] DofType " + rDof + " has no enum equivalent or is not implemented.");
}


}
}
#endif /* NODEENUM_H_ */
