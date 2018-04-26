#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Triangle.h"

namespace NuTo
{
class InterpolationTriangleQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTriangleQuadratic>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 6, 1> ShapeFunctionsTriangleOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 6, 2> DerivativeShapeFunctionsTriangleOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTriangleOrder2(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsTriangleOrder2(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTriangleOrder2(nodeId);
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
    Triangle mShape;
};
} /* NuTo */
