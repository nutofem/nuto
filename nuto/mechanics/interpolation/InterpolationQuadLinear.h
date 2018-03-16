#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/math/shapes/Quadrilateral.h"

namespace NuTo
{
class InterpolationQuadLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationQuadLinear>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions2D::ShapeFunctionsQuadOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder1(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return ShapeFunctions2D::NodeCoordinatesQuadOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 4;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Quadrilateral mShape;
};
} /* NuTo */
