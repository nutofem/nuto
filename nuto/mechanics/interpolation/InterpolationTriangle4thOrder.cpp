#include "nuto/base/Exception.h"
#include "InterpolationTriangle4thOrder.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTriangle4thOrder::Clone() const
{
    return std::make_unique<InterpolationTriangle4thOrder>(*this);
}

ShapeFunctions InterpolationTriangle4thOrder::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctionsTriangleOrder4(naturalIpCoords);
}

DerivativeShapeFunctionsNatural
InterpolationTriangle4thOrder::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctionsTriangleOrder4(naturalIpCoords);
}

NaturalCoords InterpolationTriangle4thOrder::GetLocalCoords(int nodeId) const
{
    return NodeCoordinatesTriangleOrder4(nodeId);
}

int InterpolationTriangle4thOrder::GetNumNodes() const
{
    return 15;
}

const Shape& InterpolationTriangle4thOrder::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 2, 1> InterpolationTriangle4thOrder::NodeCoordinatesTriangleOrder4(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector2d(.00, .00);
    case 1:
        return Eigen::Vector2d(.25, .00);
    case 2:
        return Eigen::Vector2d(.50, .00);
    case 3:
        return Eigen::Vector2d(.75, .00);
    case 4:
        return Eigen::Vector2d(1.00, .00);
    case 5:
        return Eigen::Vector2d(.00, .25);
    case 6:
        return Eigen::Vector2d(.25, .25);
    case 7:
        return Eigen::Vector2d(.50, .25);
    case 8:
        return Eigen::Vector2d(.75, .25);
    case 9:
        return Eigen::Vector2d(.00, .50);
    case 10:
        return Eigen::Vector2d(.25, .50);
    case 11:
        return Eigen::Vector2d(.50, .50);
    case 12:
        return Eigen::Vector2d(.00, .75);
    case 13:
        return Eigen::Vector2d(.25, .75);
    case 14:
        return Eigen::Vector2d(.00, 1.00);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..14)");
    }
}

Eigen::Matrix<double, 15, 1>
InterpolationTriangle4thOrder::ShapeFunctionsTriangleOrder4(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 15, 1> shapeFunctions;
    double r(rCoordinates(0));
    double s(rCoordinates(1));

    shapeFunctions[0] = +1.0 - 8.33333333333 * r - 8.33333333333 * s + 23.3333333333 * r * r + 46.6666666667 * r * s +
                        23.3333333333 * s * s - 26.6666666667 * r * r * r - 80.0 * r * r * s - 80.0 * r * s * s -
                        26.6666666667 * s * s * s + 10.6666666667 * s * s * s * s + 42.6666666667 * r * s * s * s +
                        64.0 * r * r * s * s + 42.6666666667 * r * r * r * s + 10.6666666667 * r * r * r * r;
    shapeFunctions[1] = +16.0 * r - 69.3333333333 * r * r - 69.3333333333 * r * s + 96.0 * r * r * r +
                        192.0 * r * r * s + 96.0 * r * s * s - 42.6666666667 * r * s * s * s - 128.0 * r * r * s * s -
                        128.0 * r * r * r * s - 42.6666666667 * r * r * r * r;
    shapeFunctions[2] = -12.0 * r + 76.0 * r * r + 28.0 * r * s - 128.0 * r * r * r - 144.0 * r * r * s -
                        16.0 * r * s * s + 64.0 * r * r * s * s + 128.0 * r * r * r * s + 64.0 * r * r * r * r;
    shapeFunctions[3] = +5.33333333333 * r - 37.3333333333 * r * r - 5.33333333333 * r * s + 74.6666666667 * r * r * r +
                        32.0 * r * r * s - 42.6666666667 * r * r * r * s - 42.6666666667 * r * r * r * r;
    shapeFunctions[4] = -1.0 * r + 7.33333333333 * r * r - 16.0 * r * r * r + 10.6666666667 * r * r * r * r;
    shapeFunctions[5] = +16.0 * s - 69.3333333333 * r * s - 69.3333333333 * s * s + 96.0 * r * r * s +
                        192.0 * r * s * s + 96.0 * s * s * s - 42.6666666667 * s * s * s * s - 128.0 * r * s * s * s -
                        128.0 * r * r * s * s - 42.6666666667 * r * r * r * s;
    shapeFunctions[6] = +96.0 * r * s - 224.0 * r * r * s - 224.0 * r * s * s + 128.0 * r * s * s * s +
                        256.0 * r * r * s * s + 128.0 * r * r * r * s;
    shapeFunctions[7] =
            -32.0 * r * s + 160.0 * r * r * s + 32.0 * r * s * s - 128.0 * r * r * s * s - 128.0 * r * r * r * s;
    shapeFunctions[8] = +5.33333333333 * r * s - 32.0 * r * r * s + 42.6666666667 * r * r * r * s;
    shapeFunctions[9] = -12.0 * s + 28.0 * r * s + 76.0 * s * s - 16.0 * r * r * s - 144.0 * r * s * s -
                        128.0 * s * s * s + 64.0 * s * s * s * s + 128.0 * r * s * s * s + 64.0 * r * r * s * s;
    shapeFunctions[10] =
            -32.0 * r * s + 32.0 * r * r * s + 160.0 * r * s * s - 128.0 * r * s * s * s - 128.0 * r * r * s * s;
    shapeFunctions[11] = +4.0 * r * s - 16.0 * r * r * s - 16.0 * r * s * s + 64.0 * r * r * s * s;
    shapeFunctions[12] = +5.33333333333 * s - 5.33333333333 * r * s - 37.3333333333 * s * s + 32.0 * r * s * s +
                         74.6666666667 * s * s * s - 42.6666666667 * s * s * s * s - 42.6666666667 * r * s * s * s;
    shapeFunctions[13] = +5.33333333333 * r * s - 32.0 * r * s * s + 42.6666666667 * r * s * s * s;
    shapeFunctions[14] = -1.0 * s + 7.33333333333 * s * s - 16.0 * s * s * s + 10.6666666667 * s * s * s * s;

    return shapeFunctions;
}

Eigen::Matrix<double, 15, 2>
InterpolationTriangle4thOrder::DerivativeShapeFunctionsTriangleOrder4(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 15, 2> derivativeShapeFunctions;
    double r(rCoordinates(0));
    double s(rCoordinates(1));

    derivativeShapeFunctions(0, 0) = -8.33333333333 + 46.6666666667 * r + 46.6666666667 * s - 80.0 * r * r -
                                     160.0 * r * s - 80.0 * s * s + 42.6666666667 * s * s * s + 128.0 * r * s * s +
                                     128.0 * r * r * s + 42.6666666667 * r * r * r;
    derivativeShapeFunctions(0, 1) = -8.33333333333 + 46.6666666667 * r + 46.6666666667 * s - 80.0 * r * r -
                                     160.0 * r * s - 80.0 * s * s + 42.6666666667 * s * s * s + 128.0 * r * s * s +
                                     128.0 * r * r * s + 42.6666666667 * r * r * r;

    derivativeShapeFunctions(1, 0) = +16.0 - 138.666666667 * r - 69.3333333333 * s + 288.0 * r * r + 384.0 * r * s +
                                     96.0 * s * s - 42.6666666667 * s * s * s - 256.0 * r * s * s - 384.0 * r * r * s -
                                     170.666666667 * r * r * r;
    derivativeShapeFunctions(1, 1) = -69.3333333333 * r + 192.0 * r * r + 192.0 * r * s - 128.0 * r * s * s -
                                     256.0 * r * r * s - 128.0 * r * r * r;

    derivativeShapeFunctions(2, 0) = -12.0 + 152.0 * r + 28.0 * s - 384.0 * r * r - 288.0 * r * s - 16.0 * s * s +
                                     128.0 * r * s * s + 384.0 * r * r * s + 256.0 * r * r * r;
    derivativeShapeFunctions(2, 1) = +28.0 * r - 144.0 * r * r - 32.0 * r * s + 128.0 * r * r * s + 128.0 * r * r * r;

    derivativeShapeFunctions(3, 0) = +5.33333333333 - 74.6666666667 * r - 5.33333333333 * s + 224.0 * r * r +
                                     64.0 * r * s - 128.0 * r * r * s - 170.666666667 * r * r * r;
    derivativeShapeFunctions(3, 1) = -5.33333333333 * r + 32.0 * r * r - 42.6666666667 * r * r * r;

    derivativeShapeFunctions(4, 0) = -1.0 + 14.6666666667 * r - 48.0 * r * r + 42.6666666667 * r * r * r;
    derivativeShapeFunctions(4, 1) = 0.;

    derivativeShapeFunctions(5, 0) = -69.3333333333 * s + 192.0 * r * s + 192.0 * s * s - 128.0 * s * s * s -
                                     256.0 * r * s * s - 128.0 * r * r * s;
    derivativeShapeFunctions(5, 1) = +16.0 - 69.3333333333 * r - 138.666666667 * s + 96.0 * r * r + 384.0 * r * s +
                                     288.0 * s * s - 170.666666667 * s * s * s - 384.0 * r * s * s - 256.0 * r * r * s -
                                     42.6666666667 * r * r * r;

    derivativeShapeFunctions(6, 0) =
            +96.0 * s - 448.0 * r * s - 224.0 * s * s + 128.0 * s * s * s + 512.0 * r * s * s + 384.0 * r * r * s;
    derivativeShapeFunctions(6, 1) =
            +96.0 * r - 224.0 * r * r - 448.0 * r * s + 384.0 * r * s * s + 512.0 * r * r * s + 128.0 * r * r * r;

    derivativeShapeFunctions(7, 0) = -32.0 * s + 320.0 * r * s + 32.0 * s * s - 256.0 * r * s * s - 384.0 * r * r * s;
    derivativeShapeFunctions(7, 1) = -32.0 * r + 160.0 * r * r + 64.0 * r * s - 256.0 * r * r * s - 128.0 * r * r * r;

    derivativeShapeFunctions(8, 0) = +5.33333333333 * s - 64.0 * r * s + 128.0 * r * r * s;
    derivativeShapeFunctions(8, 1) = +5.33333333333 * r - 32.0 * r * r + 42.6666666667 * r * r * r;

    derivativeShapeFunctions(9, 0) = +28.0 * s - 32.0 * r * s - 144.0 * s * s + 128.0 * s * s * s + 128.0 * r * s * s;
    derivativeShapeFunctions(9, 1) = -12.0 + 28.0 * r + 152.0 * s - 16.0 * r * r - 288.0 * r * s - 384.0 * s * s +
                                     256.0 * s * s * s + 384.0 * r * s * s + 128.0 * r * r * s;

    derivativeShapeFunctions(10, 0) = -32.0 * s + 64.0 * r * s + 160.0 * s * s - 128.0 * s * s * s - 256.0 * r * s * s;
    derivativeShapeFunctions(10, 1) = -32.0 * r + 32.0 * r * r + 320.0 * r * s - 384.0 * r * s * s - 256.0 * r * r * s;

    derivativeShapeFunctions(11, 0) = +4.0 * s - 32.0 * r * s - 16.0 * s * s + 128.0 * r * s * s;
    derivativeShapeFunctions(11, 1) = +4.0 * r - 16.0 * r * r - 32.0 * r * s + 128.0 * r * r * s;

    derivativeShapeFunctions(12, 0) = -5.33333333333 * s + 32.0 * s * s - 42.6666666667 * s * s * s;
    derivativeShapeFunctions(12, 1) = +5.33333333333 - 5.33333333333 * r - 74.6666666667 * s + 64.0 * r * s +
                                      224.0 * s * s - 170.666666667 * s * s * s - 128.0 * r * s * s;

    derivativeShapeFunctions(13, 0) = +5.33333333333 * s - 32.0 * s * s + 42.6666666667 * s * s * s;
    derivativeShapeFunctions(13, 1) = +5.33333333333 * r - 64.0 * r * s + 128.0 * r * s * s;

    derivativeShapeFunctions(14, 0) = 0.;
    derivativeShapeFunctions(14, 1) = -1.0 + 14.6666666667 * s - 48.0 * s * s + 42.6666666667 * s * s * s;

    return derivativeShapeFunctions;
}
