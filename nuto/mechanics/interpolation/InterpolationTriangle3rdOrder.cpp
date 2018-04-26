#include "nuto/base/Exception.h"
#include "InterpolationTriangle3rdOrder.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTriangle3rdOrder::Clone() const
{
    return std::make_unique<InterpolationTriangle3rdOrder>(*this);
}

ShapeFunctions InterpolationTriangle3rdOrder::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctionsTriangleOrder3(naturalIpCoords);
}

DerivativeShapeFunctionsNatural
InterpolationTriangle3rdOrder::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctionsTriangleOrder3(naturalIpCoords);
}

NaturalCoords InterpolationTriangle3rdOrder::GetLocalCoords(int nodeId) const
{
    return NodeCoordinatesTriangleOrder3(nodeId);
}

int InterpolationTriangle3rdOrder::GetNumNodes() const
{
    return 10;
}

const Shape& InterpolationTriangle3rdOrder::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 2, 1> InterpolationTriangle3rdOrder::NodeCoordinatesTriangleOrder3(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector2d(0.0, 0.0);
    case 1:
        return Eigen::Vector2d(1. / 3., 0.0);
    case 2:
        return Eigen::Vector2d(2. / 3., 0.0);
    case 3:
        return Eigen::Vector2d(1.0, 0.0);
    case 4:
        return Eigen::Vector2d(0.0, 1. / 3.);
    case 5:
        return Eigen::Vector2d(1. / 3., 1. / 3.);
    case 6:
        return Eigen::Vector2d(2. / 3., 1. / 3.);
    case 7:
        return Eigen::Vector2d(0., 2. / 3.);
    case 8:
        return Eigen::Vector2d(1. / 3., 2. / 3.);
    case 9:
        return Eigen::Vector2d(0., 1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..9)");
    }
}

Eigen::Matrix<double, 10, 1>
InterpolationTriangle3rdOrder::ShapeFunctionsTriangleOrder3(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 10, 1> shapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];

    shapeFunctions[0] = +1.0 - 5.5 * r - 5.5 * s + 9.0 * r * r + 18.0 * r * s + 9.0 * s * s - 4.5 * r * r * r -
                        13.5 * r * r * s - 13.5 * r * s * s - 4.5 * s * s * s;
    shapeFunctions[1] = +9.0 * r - 22.5 * r * r - 22.5 * r * s + 13.5 * r * r * r + 27.0 * r * r * s + 13.5 * r * s * s;
    shapeFunctions[2] = -4.5 * r + 18.0 * r * r + 4.5 * r * s - 13.5 * r * r * r - 13.5 * r * r * s;
    shapeFunctions[3] = +1.0 * r - 4.5 * r * r + 4.5 * r * r * r;
    shapeFunctions[4] = +9.0 * s - 22.5 * r * s - 22.5 * s * s + 13.5 * r * r * s + 27.0 * r * s * s + 13.5 * s * s * s;
    shapeFunctions[5] = +27.0 * r * s - 27.0 * r * r * s - 27.0 * r * s * s;
    shapeFunctions[6] = -4.5 * r * s + 13.5 * r * r * s;
    shapeFunctions[7] = -4.5 * s + 4.5 * r * s + 18.0 * s * s - 13.5 * r * s * s - 13.5 * s * s * s;
    shapeFunctions[8] = -4.5 * r * s + 13.5 * r * s * s;
    shapeFunctions[9] = +1.0 * s - 4.5 * s * s + 4.5 * s * s * s;

    return shapeFunctions;
}

Eigen::Matrix<double, 10, 2>
InterpolationTriangle3rdOrder::DerivativeShapeFunctionsTriangleOrder3(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 10, 2> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];

    derivativeShapeFunctions(0, 0) = -5.5 + 18.0 * r + 18.0 * s - 13.5 * r * r - 27.0 * r * s - 13.5 * s * s;
    derivativeShapeFunctions(0, 1) = -5.5 + 18.0 * r + 18.0 * s - 13.5 * r * r - 27.0 * r * s - 13.5 * s * s;

    derivativeShapeFunctions(1, 0) = +9.0 - 45.0 * r - 22.5 * s + 40.5 * r * r + 54.0 * r * s + 13.5 * s * s;
    derivativeShapeFunctions(1, 1) = -22.5 * r + 27.0 * r * r + 27.0 * r * s;

    derivativeShapeFunctions(2, 0) = -4.5 + 36.0 * r + 4.5 * s - 40.5 * r * r - 27.0 * r * s;
    derivativeShapeFunctions(2, 1) = +4.5 * r - 13.5 * r * r;

    derivativeShapeFunctions(3, 0) = +1.0 - 9.0 * r + 13.5 * r * r;
    derivativeShapeFunctions(3, 1) = 0.;

    derivativeShapeFunctions(4, 0) = -22.5 * s + 27.0 * r * s + 27.0 * s * s;
    derivativeShapeFunctions(4, 1) = +9.0 - 22.5 * r - 45.0 * s + 13.5 * r * r + 54.0 * r * s + 40.5 * s * s;

    derivativeShapeFunctions(5, 0) = +27.0 * s - 54.0 * r * s - 27.0 * s * s;
    derivativeShapeFunctions(5, 1) = +27.0 * r - 27.0 * r * r - 54.0 * r * s;

    derivativeShapeFunctions(6, 0) = -4.5 * s + 27.0 * r * s;
    derivativeShapeFunctions(6, 1) = -4.5 * r + 13.5 * r * r;

    derivativeShapeFunctions(7, 0) = +4.5 * s - 13.5 * s * s;
    derivativeShapeFunctions(7, 1) = -4.5 + 4.5 * r + 36.0 * s - 27.0 * r * s - 40.5 * s * s;

    derivativeShapeFunctions(8, 0) = -4.5 * s + 13.5 * s * s;
    derivativeShapeFunctions(8, 1) = -4.5 * r + 27.0 * r * s;

    derivativeShapeFunctions(9, 0) = 0.;
    derivativeShapeFunctions(9, 1) = +1.0 - 9.0 * s + 13.5 * s * s;

    return derivativeShapeFunctions;
}
