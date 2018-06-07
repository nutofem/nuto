#include "nuto/base/Exception.h"
#include "InterpolationTruss4thOrder.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTruss4thOrder::Clone() const
{
    return std::make_unique<InterpolationTruss4thOrder>(*this);
}

Eigen::VectorXd InterpolationTruss4thOrder::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd InterpolationTruss4thOrder::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationTruss4thOrder::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationTruss4thOrder::GetNumNodes() const
{
    return 5;
}

const Shape& InterpolationTruss4thOrder::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 1, 1> InterpolationTruss4thOrder::LocalCoords(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1:
        return Eigen::Matrix<double, 1, 1>::Constant(-0.5);
    case 2:
        return Eigen::Matrix<double, 1, 1>::Constant(0.);
    case 3:
        return Eigen::Matrix<double, 1, 1>::Constant(0.5);
    case 4:
        return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..4)");
    }
}

Eigen::Matrix<double, 5, 1> InterpolationTruss4thOrder::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 5, 1> shapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r2 * r;
    double r4 = r3 * r;
    shapeFunctions[0] = +0.166666666667 * r - 0.166666666667 * r2 - 0.666666666667 * r3 + 0.666666666667 * r4;
    shapeFunctions[1] = -1.33333333333 * r + 2.66666666667 * r2 + 1.33333333333 * r3 - 2.66666666667 * r4;
    shapeFunctions[2] = +1.0 - 5.0 * r2 + 4.0 * r4;
    shapeFunctions[3] = +1.33333333333 * r + 2.66666666667 * r2 - 1.33333333333 * r3 - 2.66666666667 * r4;
    shapeFunctions[4] = -0.166666666667 * r - 0.166666666667 * r2 + 0.666666666667 * r3 + 0.666666666667 * r4;
    return shapeFunctions;
}

Eigen::Matrix<double, 5, 1> InterpolationTruss4thOrder::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 5, 1> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r2 * r;
    derivativeShapeFunctions[0] = +0.166666666667 - 0.333333333333 * r - 2.0 * r2 + 2.66666666667 * r3;
    derivativeShapeFunctions[1] = -1.33333333333 + 5.33333333333 * r + 4.0 * r2 - 10.6666666667 * r3;
    derivativeShapeFunctions[2] = -10.0 * r + 16.0 * r3;
    derivativeShapeFunctions[3] = +1.33333333333 + 5.33333333333 * r - 4.0 * r2 - 10.6666666667 * r3;
    derivativeShapeFunctions[4] = -0.166666666667 - 0.333333333333 * r + 2.0 * r2 + 2.66666666667 * r3;
    return derivativeShapeFunctions;
}

std::vector<int> InterpolationTruss4thOrder::EdgeNodeIds(int /* edgeIndex */) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

std::unique_ptr<InterpolationSimple> InterpolationTruss4thOrder::EdgeInterpolation(int /* edgeIndex*/) const
{
    return std::make_unique<InterpolationTruss4thOrder>();
}

std::vector<int> InterpolationTruss4thOrder::FaceNodeIds(int /* faceIndex */) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}

std::unique_ptr<InterpolationSimple> InterpolationTruss4thOrder::FaceInterpolation(int /* faceIndex*/) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented");
}
