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

enum eAttributes
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

static inline std::map<eAttributes, std::string> GetAttributeMap()
{
    std::map<eAttributes, std::string> attributeMap;
    attributeMap[eAttributes::COORDINATES]             = "COORDINATES";
    attributeMap[eAttributes::ROTATIONS]               = "ROTATIONS";
    attributeMap[eAttributes::TEMPERATURES]            = "TEMPERATURES";
    attributeMap[eAttributes::DISPLACEMENTS]           = "DISPLACEMENTS";
    attributeMap[eAttributes::FINESCALEDISPLACEMENTS]  = "FINESCALEDISPLACEMENTS";
    attributeMap[eAttributes::NONLOCALDATA]            = "NONLOCALDATA";
    attributeMap[eAttributes::NONLOCALEQPLASTICSTRAIN] = "NONLOCALEQPLASTICSTRAIN";
    attributeMap[eAttributes::NONLOCALTOTALSTRAIN]     = "NONLOCALTOTALSTRAIN";
    attributeMap[eAttributes::NONLOCALEQSTRAIN]        = "NONLOCALEQSTRAIN";
    attributeMap[eAttributes::WATERVOLUMEFRACTION]     = "WATERVOLUMEFRACTION";
    attributeMap[eAttributes::RELATIVEHUMIDITY]        = "RELATIVEHUMIDITY";
    return attributeMap;
}


static inline std::string AttributeToString(const eAttributes& rAttribute)
{
    try
    {
        return GetAttributeMap().find(rAttribute)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException("[NuTo::Node::AttributeToString] Enum undefined or not implemented.");
    }
}

static inline eAttributes AttributeToEnum(const std::string& rAttribute)
{
    std::string uppercase = boost::to_upper_copy(rAttribute);

    for(auto entry : GetAttributeMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Node::AttributeToEnum] DofType " + rAttribute + " has no enum equivalent or is not implemented.");
}


}
}
#endif /* NODEENUM_H_ */
