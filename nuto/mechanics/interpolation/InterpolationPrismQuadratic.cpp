#include "nuto/base/Exception.h"
#include "InterpolationPrismQuadratic.h"
#include "InterpolationTriangleQuadratic.h"
#include "InterpolationTrussQuadratic.h"

using namespace NuTo;

Eigen::Matrix<double, 3, 1> InterpolationPrismQuadratic::NodeCoordinatesPrismOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0:
        return Eigen::Vector3d(0., 0., -1.);
    case 1:
        return Eigen::Vector3d(1., 0., -1.);
    case 2:
        return Eigen::Vector3d(0., 1., -1.);
    case 3:
        return Eigen::Vector3d(0., 0., 1.);
    case 4:
        return Eigen::Vector3d(1., 0., 1.);
    case 5:
        return Eigen::Vector3d(0., 1., 1.);

    case 6:
        return Eigen::Vector3d(.5, 0., -1.);
    case 7:
        return Eigen::Vector3d(0., .5, -1.);
    case 8:
        return Eigen::Vector3d(0., 0., 0.);
    case 9:
        return Eigen::Vector3d(.5, .5, -1.);
    case 10:
        return Eigen::Vector3d(1., 0., 0.);
    case 11:
        return Eigen::Vector3d(0, 1., 0.);

    case 12:
        return Eigen::Vector3d(.5, 0., 1.);
    case 13:
        return Eigen::Vector3d(0., .5, 1.);
    case 14:
        return Eigen::Vector3d(.5, .5, 1.);

    case 15:
        return Eigen::Vector3d(.5, 0., 0.);
    case 16:
        return Eigen::Vector3d(0., .5, 0.);
    case 17:
        return Eigen::Vector3d(.5, .5, 0.);

    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..14)");
    }
}

Eigen::Matrix<double, 18, 1> InterpolationPrismQuadratic::ShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 18, 1> shapeFunctions;

    auto NTriangle2 = InterpolationTriangleQuadratic::ShapeFunctionsTriangleOrder2(
            Eigen::Vector2d(rCoordinates[0], rCoordinates[1]));
    auto NTruss2 = InterpolationTrussQuadratic::ShapeFunctionsTrussOrder2(
            Eigen::Matrix<double, 1, 1>::Constant(rCoordinates[2]));

    shapeFunctions[0] = NTriangle2[0] * NTruss2[0];
    shapeFunctions[1] = NTriangle2[1] * NTruss2[0];
    shapeFunctions[2] = NTriangle2[2] * NTruss2[0];
    shapeFunctions[3] = NTriangle2[0] * NTruss2[2];
    shapeFunctions[4] = NTriangle2[1] * NTruss2[2];
    shapeFunctions[5] = NTriangle2[2] * NTruss2[2];
    shapeFunctions[6] = NTriangle2[3] * NTruss2[0];
    shapeFunctions[7] = NTriangle2[5] * NTruss2[0];
    shapeFunctions[8] = NTriangle2[0] * NTruss2[1];
    shapeFunctions[9] = NTriangle2[4] * NTruss2[0];
    shapeFunctions[10] = NTriangle2[1] * NTruss2[1];
    shapeFunctions[11] = NTriangle2[2] * NTruss2[1];
    shapeFunctions[12] = NTriangle2[3] * NTruss2[2];
    shapeFunctions[13] = NTriangle2[5] * NTruss2[2];
    shapeFunctions[14] = NTriangle2[4] * NTruss2[2];
    shapeFunctions[15] = NTriangle2[3] * NTruss2[1];
    shapeFunctions[16] = NTriangle2[5] * NTruss2[1];
    shapeFunctions[17] = NTriangle2[4] * NTruss2[1];

    return shapeFunctions;
}

Eigen::Matrix<double, 18, 3>
InterpolationPrismQuadratic::DerivativeShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates)
{

    auto NTriangle2 = InterpolationTriangleQuadratic::ShapeFunctionsTriangleOrder2(
            Eigen::Vector2d(rCoordinates[0], rCoordinates[1]));
    auto dNTriangle2 = InterpolationTriangleQuadratic::DerivativeShapeFunctionsTriangleOrder2(
            Eigen::Vector2d(rCoordinates[0], rCoordinates[1]));

    auto NTruss2 = InterpolationTrussQuadratic::ShapeFunctionsTrussOrder2(
            Eigen::Matrix<double, 1, 1>::Constant(rCoordinates[2]));
    auto dNTruss2 = InterpolationTrussQuadratic::DerivativeShapeFunctionsTrussOrder2(
            Eigen::Matrix<double, 1, 1>::Constant(rCoordinates[2]));

    Eigen::Matrix<double, 18, 3> derivativeShapeFunctions;

    derivativeShapeFunctions(0, 0) = dNTriangle2(0, 0) * NTruss2[0];
    derivativeShapeFunctions(0, 1) = dNTriangle2(0, 1) * NTruss2[0];
    derivativeShapeFunctions(0, 2) = NTriangle2[0] * dNTruss2[0];

    derivativeShapeFunctions(1, 0) = dNTriangle2(1, 0) * NTruss2[0];
    derivativeShapeFunctions(1, 1) = dNTriangle2(1, 1) * NTruss2[0];
    derivativeShapeFunctions(1, 2) = NTriangle2[1] * dNTruss2[0];

    derivativeShapeFunctions(2, 0) = dNTriangle2(2, 0) * NTruss2[0];
    derivativeShapeFunctions(2, 1) = dNTriangle2(2, 1) * NTruss2[0];
    derivativeShapeFunctions(2, 2) = NTriangle2[2] * dNTruss2[0];

    derivativeShapeFunctions(3, 0) = dNTriangle2(0, 0) * NTruss2[2];
    derivativeShapeFunctions(3, 1) = dNTriangle2(0, 1) * NTruss2[2];
    derivativeShapeFunctions(3, 2) = NTriangle2[0] * dNTruss2[2];

    derivativeShapeFunctions(4, 0) = dNTriangle2(1, 0) * NTruss2[2];
    derivativeShapeFunctions(4, 1) = dNTriangle2(1, 1) * NTruss2[2];
    derivativeShapeFunctions(4, 2) = NTriangle2[1] * dNTruss2[2];

    derivativeShapeFunctions(5, 0) = dNTriangle2(2, 0) * NTruss2[2];
    derivativeShapeFunctions(5, 1) = dNTriangle2(2, 1) * NTruss2[2];
    derivativeShapeFunctions(5, 2) = NTriangle2[2] * dNTruss2[2];

    derivativeShapeFunctions(6, 0) = dNTriangle2(3, 0) * NTruss2[0];
    derivativeShapeFunctions(6, 1) = dNTriangle2(3, 1) * NTruss2[0];
    derivativeShapeFunctions(6, 2) = NTriangle2[3] * dNTruss2[0];

    derivativeShapeFunctions(7, 0) = dNTriangle2(5, 0) * NTruss2[0];
    derivativeShapeFunctions(7, 1) = dNTriangle2(5, 1) * NTruss2[0];
    derivativeShapeFunctions(7, 2) = NTriangle2[5] * dNTruss2[0];

    derivativeShapeFunctions(8, 0) = dNTriangle2(0, 0) * NTruss2[1];
    derivativeShapeFunctions(8, 1) = dNTriangle2(0, 1) * NTruss2[1];
    derivativeShapeFunctions(8, 2) = NTriangle2[0] * dNTruss2[1];

    derivativeShapeFunctions(9, 0) = dNTriangle2(4, 0) * NTruss2[0];
    derivativeShapeFunctions(9, 1) = dNTriangle2(4, 1) * NTruss2[0];
    derivativeShapeFunctions(9, 2) = NTriangle2[4] * dNTruss2[0];

    derivativeShapeFunctions(10, 0) = dNTriangle2(1, 0) * NTruss2[1];
    derivativeShapeFunctions(10, 1) = dNTriangle2(1, 1) * NTruss2[1];
    derivativeShapeFunctions(10, 2) = NTriangle2[1] * dNTruss2[1];

    derivativeShapeFunctions(11, 0) = dNTriangle2(2, 0) * NTruss2[1];
    derivativeShapeFunctions(11, 1) = dNTriangle2(2, 1) * NTruss2[1];
    derivativeShapeFunctions(11, 2) = NTriangle2[2] * dNTruss2[1];

    derivativeShapeFunctions(12, 0) = dNTriangle2(3, 0) * NTruss2[2];
    derivativeShapeFunctions(12, 1) = dNTriangle2(3, 1) * NTruss2[2];
    derivativeShapeFunctions(12, 2) = NTriangle2[3] * dNTruss2[2];

    derivativeShapeFunctions(13, 0) = dNTriangle2(5, 0) * NTruss2[2];
    derivativeShapeFunctions(13, 1) = dNTriangle2(5, 1) * NTruss2[2];
    derivativeShapeFunctions(13, 2) = NTriangle2[5] * dNTruss2[2];

    derivativeShapeFunctions(14, 0) = dNTriangle2(4, 0) * NTruss2[2];
    derivativeShapeFunctions(14, 1) = dNTriangle2(4, 1) * NTruss2[2];
    derivativeShapeFunctions(14, 2) = NTriangle2[4] * dNTruss2[2];

    derivativeShapeFunctions(15, 0) = dNTriangle2(3, 0) * NTruss2[1];
    derivativeShapeFunctions(15, 1) = dNTriangle2(3, 1) * NTruss2[1];
    derivativeShapeFunctions(15, 2) = NTriangle2[3] * dNTruss2[1];

    derivativeShapeFunctions(16, 0) = dNTriangle2(5, 0) * NTruss2[1];
    derivativeShapeFunctions(16, 1) = dNTriangle2(5, 1) * NTruss2[1];
    derivativeShapeFunctions(16, 2) = NTriangle2[5] * dNTruss2[1];

    derivativeShapeFunctions(17, 0) = dNTriangle2(4, 0) * NTruss2[1];
    derivativeShapeFunctions(17, 1) = dNTriangle2(4, 1) * NTruss2[1];
    derivativeShapeFunctions(17, 2) = NTriangle2[4] * dNTruss2[1];

    return derivativeShapeFunctions;
}
