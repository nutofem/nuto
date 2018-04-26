#include "nuto/base/Exception.h"
#include "InterpolationPrismLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationPrismLinear::Clone() const
{
    return std::make_unique<InterpolationPrismLinear>(*this);
}

ShapeFunctions InterpolationPrismLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctionsPrismOrder1(naturalIpCoords);
}

DerivativeShapeFunctionsNatural
InterpolationPrismLinear::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctionsPrismOrder1(naturalIpCoords);
}

NaturalCoords InterpolationPrismLinear::GetLocalCoords(int nodeId) const
{
    return NodeCoordinatesPrismOrder1(nodeId);
}

int InterpolationPrismLinear::GetNumNodes() const
{
    return 6;
}

const Shape& InterpolationPrismLinear::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 3, 1> InterpolationPrismLinear::NodeCoordinatesPrismOrder1(int rNodeIndex)
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
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..5)");
    }
}

Eigen::Matrix<double, 6, 1> InterpolationPrismLinear::ShapeFunctionsPrismOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 6, 1> shapeFunctions;

    double triangle0 = 1. - rCoordinates[0] - rCoordinates[1];
    double triangle1 = rCoordinates[0];
    double triangle2 = rCoordinates[1];

    shapeFunctions[0] = 0.5 * triangle0 * (1. - rCoordinates[2]);
    shapeFunctions[1] = 0.5 * triangle1 * (1. - rCoordinates[2]);
    shapeFunctions[2] = 0.5 * triangle2 * (1. - rCoordinates[2]);
    shapeFunctions[3] = 0.5 * triangle0 * (1. + rCoordinates[2]);
    shapeFunctions[4] = 0.5 * triangle1 * (1. + rCoordinates[2]);
    shapeFunctions[5] = 0.5 * triangle2 * (1. + rCoordinates[2]);
    return shapeFunctions;
}

Eigen::Matrix<double, 6, 3>
InterpolationPrismLinear::DerivativeShapeFunctionsPrismOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 6, 3> derivativeShapeFunctions;

    derivativeShapeFunctions(0, 0) = -0.5 * (1. - rCoordinates[2]);
    derivativeShapeFunctions(0, 1) = -0.5 * (1. - rCoordinates[2]);
    derivativeShapeFunctions(0, 2) = -0.5 * (1. - rCoordinates[0] - rCoordinates[1]);

    derivativeShapeFunctions(1, 0) = 0.5 * (1. - rCoordinates[2]);
    derivativeShapeFunctions(1, 1) = 0.;
    derivativeShapeFunctions(1, 2) = -0.5 * rCoordinates[0];

    derivativeShapeFunctions(2, 0) = 0.;
    derivativeShapeFunctions(2, 1) = 0.5 * (1. - rCoordinates[2]);
    derivativeShapeFunctions(2, 2) = -0.5 * rCoordinates[1];

    derivativeShapeFunctions(3, 0) = -0.5 * (1. + rCoordinates[2]);
    derivativeShapeFunctions(3, 1) = -0.5 * (1. + rCoordinates[2]);
    derivativeShapeFunctions(3, 2) = 0.5 * (1. - rCoordinates[0] - rCoordinates[1]);

    derivativeShapeFunctions(4, 0) = 0.5 * (1. + rCoordinates[2]);
    derivativeShapeFunctions(4, 1) = 0.;
    derivativeShapeFunctions(4, 2) = 0.5 * rCoordinates[0];

    derivativeShapeFunctions(5, 0) = 0.;
    derivativeShapeFunctions(5, 1) = 0.5 * (1. + rCoordinates[2]);
    derivativeShapeFunctions(5, 2) = 0.5 * rCoordinates[1];

    return derivativeShapeFunctions;
}
