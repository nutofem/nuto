#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/math/shapes/Prism.h"

namespace NuTo
{
class InterpolationPrismLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationPrismLinear>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions3D::ShapeFunctionsPrismOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions3D::DerivativeShapeFunctionsPrismOrder1(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return ShapeFunctions3D::NodeCoordinatesPrismOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 6;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Prism mShape;
};
} /* NuTo */
