#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/elements/ElementShapeFunctions.h"

namespace NuTo
{
class InterpolationTetrahedronQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTetrahedronQuadratic>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions3D::ShapeFunctionsTetrahedronOrder2(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& coords) const override
    {
        return ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder2(coords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return ShapeFunctions3D::NodeCoordinatesTetrahedronOrder2(nodeId);
    }

    int GetNumNodes() const override
    {
        return 10;
    }
};
} /* NuTo */
