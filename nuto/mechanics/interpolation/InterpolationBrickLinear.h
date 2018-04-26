#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Hexahedron.h"

namespace NuTo
{
class InterpolationBrickLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationBrickLinear>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 8, 1> ShapeFunctionsBrickOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 8, 3> DerivativeShapeFunctionsBrickOrder1(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsBrickOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsBrickOrder1(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesBrickOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 8;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Hexahedron mShape;
};
} /* NuTo */
