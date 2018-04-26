#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Triangle.h"

namespace NuTo
{
class InterpolationTriangle3rdOrder : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTriangle3rdOrder>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder3(int rNodeIndex);

    static Eigen::Matrix<double, 10, 1> ShapeFunctionsTriangleOrder3(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 10, 2> DerivativeShapeFunctionsTriangleOrder3(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTriangleOrder3(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsTriangleOrder3(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTriangleOrder3(nodeId);
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
    Triangle mShape;
};
} /* NuTo */
