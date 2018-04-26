#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Tetrahedron.h"

namespace NuTo
{
class InterpolationTetrahedronQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTetrahedronQuadratic>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesTetrahedronOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 10, 1> ShapeFunctionsTetrahedronOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 10, 3> DerivativeShapeFunctionsTetrahedronOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTetrahedronOrder2(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& coords) const override
    {
        return DerivativeShapeFunctionsTetrahedronOrder2(coords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTetrahedronOrder2(nodeId);
    }

    int GetNumNodes() const override
    {
        return 10;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Tetrahedron mShape;
};
} /* NuTo */
