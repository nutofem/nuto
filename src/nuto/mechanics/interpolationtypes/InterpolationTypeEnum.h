/*
 * InterpolationTypeEnum.h
 *
 *  Created on: 30 Mar 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATIONTYPEENUM_H_
#define INTERPOLATIONTYPEENUM_H_

#include <map>
#include <boost/algorithm/string.hpp>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{

namespace Interpolation
{

enum eShapeType
{
    SPRING,
    TRUSS1D,
    TRUSSXD,
    TRIANGLE2D,
    QUAD2D,
    TETRAHEDRON3D,
    BRICK3D,
    INTERFACE,
    IGA1D,
    IGA2D
};
static inline std::map<eShapeType, std::string> GetShapeTypeMap()
{
    std::map<eShapeType, std::string> shapeTypeMap;
    shapeTypeMap[eShapeType::SPRING]        = "SPRING";
    shapeTypeMap[eShapeType::TRUSS1D]       = "TRUSS1D";
    shapeTypeMap[eShapeType::TRUSSXD]       = "TRUSSXD";
    shapeTypeMap[eShapeType::TRIANGLE2D]    = "TRIANGLE2D";
    shapeTypeMap[eShapeType::QUAD2D]        = "QUAD2D";
    shapeTypeMap[eShapeType::TETRAHEDRON3D] = "TETRAHEDRON3D";
    shapeTypeMap[eShapeType::BRICK3D]       = "BRICK3D";
    shapeTypeMap[eShapeType::INTERFACE]     = "INTERFACE";
    shapeTypeMap[eShapeType::IGA1D]         = "IGA1D";
    shapeTypeMap[eShapeType::IGA2D]         = "IGA2D";
    return shapeTypeMap;
}

enum eTypeOrder
{
    EQUIDISTANT1,
    EQUIDISTANT2,
    EQUIDISTANT3,
    EQUIDISTANT4,
    LOBATTO2,
    LOBATTO3,
    LOBATTO4,
    SPLINE
};
static inline std::map<eTypeOrder, std::string> GetTypeOrderMap()
{
    std::map<eTypeOrder, std::string> typeOrderMap;
    typeOrderMap[eTypeOrder::EQUIDISTANT1]  = "EQUIDISTANT1";
    typeOrderMap[eTypeOrder::EQUIDISTANT2]  = "EQUIDISTANT2";
    typeOrderMap[eTypeOrder::EQUIDISTANT3]  = "EQUIDISTANT3";
    typeOrderMap[eTypeOrder::EQUIDISTANT4]  = "EQUIDISTANT4";
    typeOrderMap[eTypeOrder::LOBATTO2]      = "LOBATTO2";
    typeOrderMap[eTypeOrder::LOBATTO3]      = "LOBATTO3";
    typeOrderMap[eTypeOrder::LOBATTO4]      = "LOBATTO4";
    return typeOrderMap;
}





static inline std::string ShapeTypeToString(const eShapeType& rShapeType)
{
    try
    {
        return GetShapeTypeMap().find(rShapeType)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException("[NuTo::Interpolation::ShapeTypeToString] Enum undefined or not implemented.");
    }
}


static inline std::string TypeOrderToString(const eTypeOrder& rTypeOrder)
{
    try
    {
        return GetTypeOrderMap().find(rTypeOrder)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException("[NuTo::Interpolation::TypeOrderToString] Enum undefined or not implemented.");
    }
}

static inline eShapeType ShapeTypeToEnum(const std::string& rShapeType)
{
    std::string uppercase = boost::to_upper_copy(rShapeType);

    for(auto entry : GetShapeTypeMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Interpolation::ShapeTypeToEnum] ShapeType " + rShapeType + " has no enum equivalent or is not implemented.");
}

static inline eTypeOrder TypeOrderToEnum(const std::string& rTypeOrder)
{
    std::string uppercase = boost::to_upper_copy(rTypeOrder);

    for(auto entry : GetTypeOrderMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Interpolation::TypeOrderToEnum] TypeOrder " + rTypeOrder + " has no enum equivalent or is not implemented.");
}

} // namespace Interpolation
} // namespace NuTo

#endif /* INTERPOLATIONTYPEENUM_H_ */
