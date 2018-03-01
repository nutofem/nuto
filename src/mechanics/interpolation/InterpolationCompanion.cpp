#include "InterpolationCompanion.h"

#include "base/Exception.h"

#include "InterpolationTrussLinear.h"
#include "InterpolationTriangleLinear.h"
#include "InterpolationTriangleQuadratic.h"
#include "InterpolationQuadLinear.h"
#include "InterpolationQuadQuadratic.h"
#include "InterpolationTetrahedronLinear.h"
#include "InterpolationBrickLinear.h"
#include "InterpolationPrismLinear.h"
#include "InterpolationPyramidLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> S_Line(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationTrussLinear>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Line interpolation of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<InterpolationSimple> S_Triangle(int order)
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


std::unique_ptr<InterpolationSimple> S_Quadrilateral(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationQuadLinear>();
    case 2:
        return std::make_unique<InterpolationQuadQuadratic>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Quad interpolation of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<InterpolationSimple> S_Tetrahedron(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationTetrahedronLinear>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Tet interpolation of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<InterpolationSimple> S_Hex(int order)
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


std::unique_ptr<InterpolationSimple> S_Prism(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationPrismLinear>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Prism interpolation of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<InterpolationSimple> S_Pyramid(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationPyramidLinear>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Prism interpolation of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<InterpolationSimple> NuTo::CreateSerendipityInterpolation(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return S_Line(order);
    case eShape::Triangle:
        return S_Triangle(order);
    case eShape::Quadrilateral:
        return S_Quadrilateral(order);
    case eShape::Tetrahedron:
        return S_Tetrahedron(order);
    case eShape::Hexahedron:
        return S_Hex(order);
    case eShape::Prism:
        return S_Prism(order);
    case eShape::Pyramid:
        return S_Pyramid(order);
    }
}
