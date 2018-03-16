#include "InterpolationCompanion.h"

#include "base/Exception.h"

#include "InterpolationTrussLinear.h"
#include "InterpolationTriangleLinear.h"
#include "InterpolationTriangleQuadratic.h"
#include "InterpolationQuadLinear.h"
#include "InterpolationQuadQuadratic.h"
#include "InterpolationTetrahedronLinear.h"
#include "InterpolationTetrahedronQuadratic.h"
#include "InterpolationBrickLinear.h"
#include "InterpolationPrismLinear.h"
#include "InterpolationPrismQuadratic.h"
#include "InterpolationPyramidLinear.h"

#include "InterpolationTrussLobatto.h"
#include "InterpolationQuadLobatto.h"
#include "InterpolationBrickLobatto.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> LagrangeLine(int order)
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


std::unique_ptr<InterpolationSimple> LagrangeTriangle(int order)
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


std::unique_ptr<InterpolationSimple> LagrangeQuadrilateral(int order)
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


std::unique_ptr<InterpolationSimple> LagrangeTetrahedron(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationTetrahedronLinear>();
    case 2:
        return std::make_unique<InterpolationTetrahedronQuadratic>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Tet interpolation of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<InterpolationSimple> LagrangeHex(int order)
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


std::unique_ptr<InterpolationSimple> LagrangePrism(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<InterpolationPrismLinear>();
    case 2:
        return std::make_unique<InterpolationPrismQuadratic>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Prism interpolation of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<InterpolationSimple> LagrangePyramid(int order)
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


std::unique_ptr<InterpolationSimple> NuTo::CreateLagrangeInterpolation(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return LagrangeLine(order);
    case eShape::Triangle:
        return LagrangeTriangle(order);
    case eShape::Quadrilateral:
        return LagrangeQuadrilateral(order);
    case eShape::Tetrahedron:
        return LagrangeTetrahedron(order);
    case eShape::Hexahedron:
        return LagrangeHex(order);
    case eShape::Prism:
        return LagrangePrism(order);
    case eShape::Pyramid:
        return LagrangePyramid(order);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Unknown shape.");
    }
}


std::unique_ptr<InterpolationSimple> LobattoLine(int order)
{
    return std::make_unique<InterpolationTrussLobatto>(order);
}


std::unique_ptr<InterpolationSimple> LobattoQuadrilateral(int order)
{
    return std::make_unique<InterpolationQuadLobatto>(order);
}


std::unique_ptr<InterpolationSimple> LobattoHexahedron(int order)
{
    return std::make_unique<InterpolationBrickLobatto>(order);
}


std::unique_ptr<InterpolationSimple> NuTo::CreateLobattoInterpolation(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return LobattoLine(order);
    case eShape::Quadrilateral:
        return LobattoQuadrilateral(order);
    case eShape::Hexahedron:
        return LobattoHexahedron(order);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Lobatto interpolations are not defined for this shape.");
    }
}
