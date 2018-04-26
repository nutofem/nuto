#include "nuto/base/Exception.h"
#include "InterpolationPyramidLinear.h"
#include <cfloat>

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationPyramidLinear::Clone() const
{
    return std::make_unique<InterpolationPyramidLinear>(*this);
}

ShapeFunctions InterpolationPyramidLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctionsPyramidOrder1(naturalIpCoords);
}

DerivativeShapeFunctionsNatural
InterpolationPyramidLinear::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctionsPyramidOrder1(naturalIpCoords);
}

NaturalCoords InterpolationPyramidLinear::GetLocalCoords(int nodeId) const
{
    return NodeCoordinatesPyramidOrder1(nodeId);
}

int InterpolationPyramidLinear::GetNumNodes() const
{
    return 5;
}

const Shape& InterpolationPyramidLinear::GetShape() const
{
    return mShape;
}


Eigen::Matrix<double, 3, 1> InterpolationPyramidLinear::NodeCoordinatesPyramidOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0:
        return Eigen::Vector3d(1., 1., 0.);
    case 1:
        return Eigen::Vector3d(1., -1., 0.);
    case 2:
        return Eigen::Vector3d(-1., -1., 0.);
    case 3:
        return Eigen::Vector3d(-1., 1., 0.);
    case 4:
        return Eigen::Vector3d(0., 0., 1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..4)");
    }
}

Eigen::Matrix<double, 5, 1> InterpolationPyramidLinear::ShapeFunctionsPyramidOrder1(const Eigen::VectorXd& rCoordinates)
{
    const double x = rCoordinates[0];
    const double y = rCoordinates[1];
    const double z = rCoordinates[2];

    double rationalTerm = x * y * z;
    if ((1. - z) != 0.)
    {
        rationalTerm /= (1. - z);
    }

    Eigen::Matrix<double, 5, 1> shapeFunctions;

    shapeFunctions[0] = 0.25 * ((1 + x) * (1 + y) - z + rationalTerm);
    shapeFunctions[1] = 0.25 * ((1 + x) * (1 - y) - z - rationalTerm);
    shapeFunctions[2] = 0.25 * ((1 - x) * (1 - y) - z + rationalTerm);
    shapeFunctions[3] = 0.25 * ((1 - x) * (1 + y) - z - rationalTerm);
    shapeFunctions[4] = z;
    return shapeFunctions;
}

Eigen::Matrix<double, 5, 3>
InterpolationPyramidLinear::DerivativeShapeFunctionsPyramidOrder1(const Eigen::VectorXd& rCoordinates)
{
    const double x = rCoordinates[0];
    const double y = rCoordinates[1];
    double z = rCoordinates[2];

    Eigen::Matrix<double, 5, 3> derivativeShapeFunctions;

    if ((1. - z) == 0.)
    {
        z = 1 - DBL_EPSILON;
    }

    derivativeShapeFunctions(0, 0) = 0.25 * (1 + y / (1 - z));
    derivativeShapeFunctions(0, 1) = 0.25 * (1 + x / (1 - z));
    derivativeShapeFunctions(0, 2) = 0.25 * (-1 + x * y / (1 - z) / (1 - z));

    derivativeShapeFunctions(1, 0) = 0.25 * (1 - y / (1 - z));
    derivativeShapeFunctions(1, 1) = 0.25 * (-1 - x / (1 - z));
    derivativeShapeFunctions(1, 2) = 0.25 * (-1 - x * y / (1 - z) / (1 - z));

    derivativeShapeFunctions(2, 0) = 0.25 * (-1 + y / (1 - z));
    derivativeShapeFunctions(2, 1) = 0.25 * (-1 + x / (1 - z));
    derivativeShapeFunctions(2, 2) = 0.25 * (-1 + x * y / (1 - z) / (1 - z));

    derivativeShapeFunctions(3, 0) = 0.25 * (-1 - y / (1 - z));
    derivativeShapeFunctions(3, 1) = 0.25 * (1 - x / (1 - z));
    derivativeShapeFunctions(3, 2) = 0.25 * (-1 - x * y / (1 - z) / (1 - z));

    derivativeShapeFunctions(4, 0) = 0.;
    derivativeShapeFunctions(4, 1) = 0.;
    derivativeShapeFunctions(4, 2) = 1.;

    return derivativeShapeFunctions;
}
