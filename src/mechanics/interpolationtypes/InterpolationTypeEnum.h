/*
 * InterpolationTypeEnum.h
 *
 *  Created on: 30 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include <map>
#include <string>

namespace NuTo
{

namespace Interpolation
{

enum class eShapeType
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

const std::map<eShapeType, std::string> GetShapeTypeMap();
std::string ShapeTypeToString(const eShapeType& rShapeType);
eShapeType ShapeTypeToEnum(const std::string& rShapeType);


enum class eTypeOrder
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

const std::map<eTypeOrder, std::string> GetTypeOrderMap();
std::string TypeOrderToString(const eTypeOrder& rTypeOrder);
eTypeOrder TypeOrderToEnum(const std::string& rTypeOrder);

} // namespace Interpolation
} // namespace NuTo

