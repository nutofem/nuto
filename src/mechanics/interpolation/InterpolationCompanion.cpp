#include "InterpolationCompanion.h"

#include "base/Exception.h"

#include "InterpolationTriangleLinear.h"
#include "InterpolationTriangleQuadratic.h"
#include "InterpolationBrickLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> CreateTriangleInterpolation(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationTriangleLinear>();
    case 2:
        return std::make_unique<InterpolationTriangleQuadratic>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Triangle interpolation of order " + std::to_string(order) + " is not defined.");
    }
}

std::unique_ptr<InterpolationSimple> CreateHexInterpolation(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationBrickLinear>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Hex interpolation of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<InterpolationSimple> NuTo::CreateInterpolation(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Triangle:
        return CreateTriangleInterpolation(order);
    case eShape::Hexahedron:
        return CreateHexInterpolation(order);
    }
}
