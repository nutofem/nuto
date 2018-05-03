#include "nuto/base/Exception.h"
#include "InterpolationTriangleLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTriangleLinear::Clone() const
{
    return std::make_unique<InterpolationTriangleLinear>(*this);
}

Eigen::VectorXd InterpolationTriangleLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd
InterpolationTriangleLinear::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationTriangleLinear::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationTriangleLinear::GetNumNodes() const
{
    return 3;
}

const Shape& InterpolationTriangleLinear::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 2, 1> InterpolationTriangleLinear::LocalCoords(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector2d(0.0, 0.0);
    case 1:
        return Eigen::Vector2d(1.0, 0.0);
    case 2:
        return Eigen::Vector2d(0.0, 1.0);
    default:
        throw Exception(__PRETTY_FUNCTION__, "node index out of range (0..2)");
    }
}

Eigen::Matrix<double, 3, 1> InterpolationTriangleLinear::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 1> shapeFunctions;
    shapeFunctions[0] = 1. - rCoordinates(0) - rCoordinates(1);
    shapeFunctions[1] = rCoordinates(0);
    shapeFunctions[2] = rCoordinates(1);
    return shapeFunctions;
}

Eigen::Matrix<double, 3, 2> InterpolationTriangleLinear::DerivativeShapeFunctions(const Eigen::VectorXd&)
{
    Eigen::Matrix<double, 3, 2> derivativeShapeFunctions;
    derivativeShapeFunctions(0, 0) = -1.0;
    derivativeShapeFunctions(0, 1) = -1.0;

    derivativeShapeFunctions(1, 0) = 1.0;
    derivativeShapeFunctions(1, 1) = 0.0;

    derivativeShapeFunctions(2, 0) = 0.0;
    derivativeShapeFunctions(2, 1) = 1.0;
    return derivativeShapeFunctions;
}
