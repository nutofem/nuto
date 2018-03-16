#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/math/shapes/Prism.h"

namespace NuTo
{
class InterpolationPrismQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationPrismQuadratic>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions3D::ShapeFunctionsPrismOrder2(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions3D::DerivativeShapeFunctionsPrismOrder2(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return ShapeFunctions3D::NodeCoordinatesPrismOrder2(nodeId);
    }

    int GetNumNodes() const override
    {
        return 18;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Prism mShape;
};
} /* NuTo */
