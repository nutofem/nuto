#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Hexahedron.h"

namespace NuTo
{
class InterpolationBrickQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationBrickQuadratic>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 20, 1> ShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 20, 3> DerivativeShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsBrickOrder2(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsBrickOrder2(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesBrickOrder2(nodeId);
    }

    int GetNumNodes() const override
    {
        return 20;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Hexahedron mShape;
};
} /* NuTo */
