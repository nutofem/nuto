#include "nuto/base/Exception.h"
#include "InterpolationBrickLinear.h"

using namespace NuTo;

Eigen::Matrix<double, 3, 1> InterpolationBrickLinear::NodeCoordinatesBrickOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0:
        return Eigen::Vector3d(-1., -1., -1.);
    case 1:
        return Eigen::Vector3d(1., -1., -1.);
    case 2:
        return Eigen::Vector3d(1., 1., -1.);
    case 3:
        return Eigen::Vector3d(-1., 1., -1.);
    case 4:
        return Eigen::Vector3d(-1., -1., 1.);
    case 5:
        return Eigen::Vector3d(1., -1., 1.);
    case 6:
        return Eigen::Vector3d(1., 1., 1.);
    case 7:
        return Eigen::Vector3d(-1., 1., 1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..7)");
    }
}

Eigen::Matrix<double, 8, 1> InterpolationBrickLinear::ShapeFunctionsBrickOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 8, 1> shapeFunctions;

    double plus_r = 1.0 + rCoordinates[0];
    double plus_s = 1.0 + rCoordinates[1];
    double plus_t = 1.0 + rCoordinates[2];

    double minus_r = 1.0 - rCoordinates[0];
    double minus_s = 1.0 - rCoordinates[1];
    double minus_t = 1.0 - rCoordinates[2];

    shapeFunctions[0] = 0.125 * minus_r * minus_s * minus_t;
    shapeFunctions[1] = 0.125 * plus_r * minus_s * minus_t;
    shapeFunctions[2] = 0.125 * plus_r * plus_s * minus_t;
    shapeFunctions[3] = 0.125 * minus_r * plus_s * minus_t;
    shapeFunctions[4] = 0.125 * minus_r * minus_s * plus_t;
    shapeFunctions[5] = 0.125 * plus_r * minus_s * plus_t;
    shapeFunctions[6] = 0.125 * plus_r * plus_s * plus_t;
    shapeFunctions[7] = 0.125 * minus_r * plus_s * plus_t;

    return shapeFunctions;
}

Eigen::Matrix<double, 8, 3>
InterpolationBrickLinear::DerivativeShapeFunctionsBrickOrder1(const Eigen::VectorXd& rCoordinates)
{
    double plus_r = 1.0 + rCoordinates[0];
    double plus_s = 1.0 + rCoordinates[1];
    double plus_t = 1.0 + rCoordinates[2];

    double minus_r = 1.0 - rCoordinates[0];
    double minus_s = 1.0 - rCoordinates[1];
    double minus_t = 1.0 - rCoordinates[2];

    Eigen::Matrix<double, 8, 3> derivativeShapeFunctions;

    derivativeShapeFunctions(0, 0) = -0.125 * minus_s * minus_t;
    derivativeShapeFunctions(0, 1) = -0.125 * minus_r * minus_t;
    derivativeShapeFunctions(0, 2) = -0.125 * minus_r * minus_s;

    derivativeShapeFunctions(1, 0) = 0.125 * minus_s * minus_t;
    derivativeShapeFunctions(1, 1) = -0.125 * plus_r * minus_t;
    derivativeShapeFunctions(1, 2) = -0.125 * plus_r * minus_s;

    derivativeShapeFunctions(2, 0) = 0.125 * plus_s * minus_t;
    derivativeShapeFunctions(2, 1) = 0.125 * plus_r * minus_t;
    derivativeShapeFunctions(2, 2) = -0.125 * plus_r * plus_s;

    derivativeShapeFunctions(3, 0) = -0.125 * plus_s * minus_t;
    derivativeShapeFunctions(3, 1) = 0.125 * minus_r * minus_t;
    derivativeShapeFunctions(3, 2) = -0.125 * minus_r * plus_s;

    derivativeShapeFunctions(4, 0) = -0.125 * minus_s * plus_t;
    derivativeShapeFunctions(4, 1) = -0.125 * minus_r * plus_t;
    derivativeShapeFunctions(4, 2) = 0.125 * minus_r * minus_s;

    derivativeShapeFunctions(5, 0) = 0.125 * minus_s * plus_t;
    derivativeShapeFunctions(5, 1) = -0.125 * plus_r * plus_t;
    derivativeShapeFunctions(5, 2) = 0.125 * plus_r * minus_s;

    derivativeShapeFunctions(6, 0) = 0.125 * plus_s * plus_t;
    derivativeShapeFunctions(6, 1) = 0.125 * plus_r * plus_t;
    derivativeShapeFunctions(6, 2) = 0.125 * plus_r * plus_s;

    derivativeShapeFunctions(7, 0) = -0.125 * plus_s * plus_t;
    derivativeShapeFunctions(7, 1) = 0.125 * minus_r * plus_t;
    derivativeShapeFunctions(7, 2) = 0.125 * minus_r * plus_s;

    return derivativeShapeFunctions;
}
