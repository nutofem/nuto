#include "InterpolationTypeEnum.h"

#include <boost/algorithm/string.hpp>
#include "mechanics/MechanicsException.h"

const std::map<NuTo::Interpolation::eShapeType, std::string> NuTo::Interpolation::GetShapeTypeMap()
{
    const std::map<eShapeType, std::string> shapeTypeMap =
       {{eShapeType::SPRING,        "SPRING"},
        {eShapeType::TRUSS1D,       "TRUSS1D"},
        {eShapeType::TRUSSXD,       "TRUSSXD"},
        {eShapeType::TRIANGLE2D,    "TRIANGLE2D"},
        {eShapeType::QUAD2D,        "QUAD2D"},
        {eShapeType::TETRAHEDRON3D, "TETRAHEDRON3D"},
        {eShapeType::BRICK3D,       "BRICK3D"},
        {eShapeType::INTERFACE,     "INTERFACE"},
        {eShapeType::IGA1D,         "IGA1D"},
        {eShapeType::IGA2D,         "IGA2D"}};
    return shapeTypeMap;
}

const std::map<NuTo::Interpolation::eTypeOrder, std::string> NuTo::Interpolation::GetTypeOrderMap()
{
    const std::map<eTypeOrder, std::string> typeOrderMap =
       {{eTypeOrder::EQUIDISTANT1,  "EQUIDISTANT1"},
        {eTypeOrder::EQUIDISTANT2,  "EQUIDISTANT2"},
        {eTypeOrder::EQUIDISTANT3,  "EQUIDISTANT3"},
        {eTypeOrder::EQUIDISTANT4,  "EQUIDISTANT4"},
        {eTypeOrder::LOBATTO2,      "LOBATTO2"},
        {eTypeOrder::LOBATTO3,      "LOBATTO3"},
        {eTypeOrder::LOBATTO4,      "LOBATTO4"},
        {eTypeOrder::SPLINE,        "SPLINE"}};
    return typeOrderMap;
}

std::string NuTo::Interpolation::ShapeTypeToString(const NuTo::Interpolation::eShapeType &rShapeType)
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

std::string NuTo::Interpolation::TypeOrderToString(const NuTo::Interpolation::eTypeOrder &rTypeOrder)
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

NuTo::Interpolation::eShapeType NuTo::Interpolation::ShapeTypeToEnum(const std::string &rShapeType)
{
    std::string uppercase = boost::to_upper_copy(rShapeType);

    for(auto entry : GetShapeTypeMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Interpolation::ShapeTypeToEnum] ShapeType " + rShapeType + " has no enum equivalent or is not implemented.");
}

NuTo::Interpolation::eTypeOrder NuTo::Interpolation::TypeOrderToEnum(const std::string &rTypeOrder)
{
    std::string uppercase = boost::to_upper_copy(rTypeOrder);

    for(auto entry : GetTypeOrderMap())
        if (entry.second == uppercase)
            return entry.first;

    throw NuTo::MechanicsException("[NuTo::Interpolation::TypeOrderToEnum] TypeOrder " + rTypeOrder + " has no enum equivalent or is not implemented.");
}
