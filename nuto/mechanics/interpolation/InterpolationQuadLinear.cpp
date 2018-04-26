#include "nuto/base/Exception.h"
#include "InterpolationQuadLinear.h"

using namespace NuTo;

Eigen::Matrix<double, 2, 1> InterpolationQuadLinear::NodeCoordinatesQuadOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector2d(-1.0, -1.0);
    case 1:
        return Eigen::Vector2d(1.0, -1.0);
    case 2:
        return Eigen::Vector2d(1.0, 1.0);
    case 3:
        return Eigen::Vector2d(-1.0, 1.0);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..3)");
    }
}

Eigen::Matrix<double, 4, 1> InterpolationQuadLinear::ShapeFunctionsQuadOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> shapeFunctions;
    shapeFunctions[0] = 0.25 * (1. - rCoordinates(0)) * (1. - rCoordinates(1));
    shapeFunctions[1] = 0.25 * (1. + rCoordinates(0)) * (1. - rCoordinates(1));
    shapeFunctions[2] = 0.25 * (1. + rCoordinates(0)) * (1. + rCoordinates(1));
    shapeFunctions[3] = 0.25 * (1. - rCoordinates(0)) * (1. + rCoordinates(1));
    return shapeFunctions;
}

Eigen::Matrix<double, 4, 2> InterpolationQuadLinear::DerivativeShapeFunctionsQuadOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 2> derivativeShapeFunctions;
    derivativeShapeFunctions(0, 0) = -0.25 * (1. - rCoordinates(1));
    derivativeShapeFunctions(0, 1) = -0.25 * (1. - rCoordinates(0));

    derivativeShapeFunctions(1, 0) = +0.25 * (1. - rCoordinates(1));
    derivativeShapeFunctions(1, 1) = -0.25 * (1. + rCoordinates(0));

    derivativeShapeFunctions(2, 0) = +0.25 * (1. + rCoordinates(1));
    derivativeShapeFunctions(2, 1) = +0.25 * (1. + rCoordinates(0));

    derivativeShapeFunctions(3, 0) = -0.25 * (1. + rCoordinates(1));
    derivativeShapeFunctions(3, 1) = +0.25 * (1. - rCoordinates(0));

    return derivativeShapeFunctions;
}
