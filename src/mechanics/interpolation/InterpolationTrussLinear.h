#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/elements/ElementShapeFunctions.h"

namespace NuTo
{
class InterpolationTrussLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTrussLinear>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions1D::ShapeFunctionsTrussOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords&) const override
    {
        return ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1();
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return ShapeFunctions1D::NodeCoordinatesTrussOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 2;
    }
};
} /* NuTo */
