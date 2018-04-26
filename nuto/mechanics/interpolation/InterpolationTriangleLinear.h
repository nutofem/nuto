#pragma once
#include <Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Triangle.h"

namespace NuTo
{
class InterpolationTriangleLinear : public InterpolationSimple
{
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 3, 1> ShapeFunctionsTriangleOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 3, 2> DerivativeShapeFunctionsTriangleOrder1(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTriangleLinear>(*this);
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTriangleOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsTriangleOrder1(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTriangleOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 3;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Triangle mShape;
};
} /* NuTo */
