#include "nuto/base/Exception.h"
#include "InterpolationTruss3rdOrder.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTruss3rdOrder::Clone() const
{
    return std::make_unique<InterpolationTruss3rdOrder>(*this);
}

Eigen::VectorXd InterpolationTruss3rdOrder::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd
InterpolationTruss3rdOrder::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationTruss3rdOrder::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationTruss3rdOrder::GetNumNodes() const
{
    return 4;
}

const Shape& InterpolationTruss3rdOrder::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 1, 1> InterpolationTruss3rdOrder::LocalCoords(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1:
        return Eigen::Matrix<double, 1, 1>::Constant(-1. / 3.);
    case 2:
        return Eigen::Matrix<double, 1, 1>::Constant(1. / 3.);
    case 3:
        return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..3)");
    }
}

Eigen::Matrix<double, 4, 1> InterpolationTruss3rdOrder::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> shapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r2 * r;
    shapeFunctions[0] = -0.0625 + 0.0625 * r + 0.5625 * r2 - 0.5625 * r3;
    shapeFunctions[1] = +0.5625 - 1.6875 * r - 0.5625 * r2 + 1.6875 * r3;
    shapeFunctions[2] = +0.5625 + 1.6875 * r - 0.5625 * r2 - 1.6875 * r3;
    shapeFunctions[3] = -0.0625 - 0.0625 * r + 0.5625 * r2 + 0.5625 * r3;
    return shapeFunctions;
}

Eigen::Matrix<double, 4, 1> InterpolationTruss3rdOrder::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    derivativeShapeFunctions[0] = +0.0625 + 1.125 * r - 1.6875 * r2;
    derivativeShapeFunctions[1] = -1.6875 - 1.125 * r + 5.0625 * r2;
    derivativeShapeFunctions[2] = +1.6875 - 1.125 * r - 5.0625 * r2;
    derivativeShapeFunctions[3] = -0.0625 + 1.125 * r + 1.6875 * r2;
    return derivativeShapeFunctions;
}
