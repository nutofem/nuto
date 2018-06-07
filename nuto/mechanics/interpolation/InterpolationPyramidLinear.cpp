#include "nuto/base/Exception.h"
#include "InterpolationPyramidLinear.h"
#include "InterpolationTriangleLinear.h"
#include "InterpolationQuadLinear.h"
#include "InterpolationTrussLinear.h"
#include <cfloat>

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationPyramidLinear::Clone() const
{
    return std::make_unique<InterpolationPyramidLinear>(*this);
}

Eigen::VectorXd InterpolationPyramidLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd InterpolationPyramidLinear::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationPyramidLinear::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationPyramidLinear::GetNumNodes() const
{
    return 5;
}

const Shape& InterpolationPyramidLinear::GetShape() const
{
    return mShape;
}


Eigen::Matrix<double, 3, 1> InterpolationPyramidLinear::LocalCoords(int rNodeIndex)
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

Eigen::Matrix<double, 5, 1> InterpolationPyramidLinear::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
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

Eigen::Matrix<double, 5, 3> InterpolationPyramidLinear::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
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

std::vector<int> InterpolationPyramidLinear::EdgeNodeIds(int edgeIndex) const
{
    switch (edgeIndex)
    {
    case 0:
        return {0, 1};
    case 1:
        return {1, 2};
    case 2:
        return {2, 3};
    case 3:
        return {3, 0};
    case 4:
        return {0, 4};
    case 5:
        return {1, 4};
    case 6:
        return {2, 4};
    case 7:
        return {3, 4};
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "edge index out of range (0..7)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationPyramidLinear::EdgeInterpolation(int /* edgeIndex*/) const
{
    return std::make_unique<InterpolationTrussLinear>();
}

std::vector<int> InterpolationPyramidLinear::FaceNodeIds(int faceIndex) const
{
    switch (faceIndex)
    {
    case 0:
        return {0, 1, 2, 3};
    case 1:
        return {0, 4, 1};
    case 2:
        return {1, 4, 2};
    case 3:
        return {2, 4, 3};
    case 4:
        return {0, 3, 4};
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "edge index out of range (0..4)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationPyramidLinear::FaceInterpolation(int faceIndex) const
{
    switch (faceIndex)
    {
    case 0:
        return std::make_unique<InterpolationQuadLinear>();
    case 1:
    case 2:
    case 3:
    case 4:
        return std::make_unique<InterpolationTriangleLinear>();
    default:
        throw Exception(__PRETTY_FUNCTION__, "Face index out of range (0-5).");
    }
}
