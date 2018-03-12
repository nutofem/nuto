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

#include "InterpolationTrussLobatto.h"
#include "InterpolationQuadLobatto.h"
#include "InterpolationBrickLobatto.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> SerendipityLine(int order)
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


std::unique_ptr<InterpolationSimple> SerendipityTriangle(int order)
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


std::unique_ptr<InterpolationSimple> SerendipityQuadrilateral(int order)
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


std::unique_ptr<InterpolationSimple> SerendipityTetrahedron(int order)
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


std::unique_ptr<InterpolationSimple> SerendipityHex(int order)
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


std::unique_ptr<InterpolationSimple> SerendipityPrism(int order)
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


std::unique_ptr<InterpolationSimple> SerendipityPyramid(int order)
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
        return SerendipityLine(order);
    case eShape::Triangle:
        return SerendipityTriangle(order);
    case eShape::Quadrilateral:
        return SerendipityQuadrilateral(order);
    case eShape::Tetrahedron:
        return SerendipityTetrahedron(order);
    case eShape::Hexahedron:
        return SerendipityHex(order);
    case eShape::Prism:
        return SerendipityPrism(order);
    case eShape::Pyramid:
        return SerendipityPyramid(order);
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
