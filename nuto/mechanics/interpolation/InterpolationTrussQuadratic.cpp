#include "nuto/base/Exception.h"
#include "InterpolationTrussQuadratic.h"

using namespace NuTo;

Eigen::Matrix<double, 1, 1> InterpolationTrussQuadratic::NodeCoordinatesTrussOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1:
        return Eigen::Matrix<double, 1, 1>::Constant(0.);
    case 2:
        return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..2)");
    }
}

Eigen::Matrix<double, 3, 1> InterpolationTrussQuadratic::ShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 1> shapeFunctions;
    shapeFunctions[0] = 0.5 * (1. - rCoordinates[0]) - 0.5 * (1. - rCoordinates[0] * rCoordinates[0]);
    shapeFunctions[1] = 1. - rCoordinates[0] * rCoordinates[0];
    shapeFunctions[2] = 0.5 * (1. + rCoordinates[0]) - 0.5 * (1. - rCoordinates[0] * rCoordinates[0]);
    return shapeFunctions;
}

Eigen::Matrix<double, 3, 1>
InterpolationTrussQuadratic::DerivativeShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 1> derivativeShapeFunctions;
    derivativeShapeFunctions[0] = -0.5 + rCoordinates[0];
    derivativeShapeFunctions[1] = -2.0 * rCoordinates[0];
    derivativeShapeFunctions[2] = 0.5 + rCoordinates[0];
    return derivativeShapeFunctions;
}
