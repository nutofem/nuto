#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/math/shapes/Tetrahedron.h"

namespace NuTo
{
class InterpolationTetrahedronLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTetrahedronLinear>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctions3D::ShapeFunctionsTetrahedronOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords&) const override
    {
        return ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder1();
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return ShapeFunctions3D::NodeCoordinatesTetrahedronOrder1(nodeId);
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
    Tetrahedron mShape;
};
} /* NuTo */
