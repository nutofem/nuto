#include "nuto/base/Exception.h"
#include "InterpolationTrussLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTrussLinear::Clone() const
{
    return std::make_unique<InterpolationTrussLinear>(*this);
}

ShapeFunctions InterpolationTrussLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctionsTrussOrder1(naturalIpCoords);
}

DerivativeShapeFunctionsNatural InterpolationTrussLinear::GetDerivativeShapeFunctions(const NaturalCoords&) const
{
    return DerivativeShapeFunctionsTrussOrder1();
}

NaturalCoords InterpolationTrussLinear::GetLocalCoords(int nodeId) const
{
    return NodeCoordinatesTrussOrder1(nodeId);
}

int InterpolationTrussLinear::GetNumNodes() const
{
    return 2;
}

const Shape& InterpolationTrussLinear::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 1, 1> InterpolationTrussLinear::NodeCoordinatesTrussOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1:
        return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..1)");
    }
}

Eigen::Matrix<double, 2, 1> InterpolationTrussLinear::ShapeFunctionsTrussOrder1(const Eigen::VectorXd& rCoordinates)
{
    return Eigen::Vector2d(0.5 * (1. - rCoordinates[0]), 0.5 * (1. + rCoordinates[0]));
}

Eigen::Matrix<double, 2, 1> InterpolationTrussLinear::DerivativeShapeFunctionsTrussOrder1()
{
    return Eigen::Vector2d(-0.5, 0.5);
}
