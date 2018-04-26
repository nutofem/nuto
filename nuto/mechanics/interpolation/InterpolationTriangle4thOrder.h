#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Triangle.h"

namespace NuTo
{
class InterpolationTriangle4thOrder : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTriangle4thOrder>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder4(int rNodeIndex);

    static Eigen::Matrix<double, 15, 1> ShapeFunctionsTriangleOrder4(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 15, 2> DerivativeShapeFunctionsTriangleOrder4(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTriangleOrder4(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsTriangleOrder4(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTriangleOrder4(nodeId);
    }

    int GetNumNodes() const override
    {
        return 15;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Triangle mShape;
};
} /* NuTo */
