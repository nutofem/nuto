#include "nuto/base/Exception.h"
#include "InterpolationPrismLinear.h"
#include "InterpolationTrussLinear.h"
#include "InterpolationTriangleLinear.h"
#include "InterpolationQuadLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationPrismLinear::Clone() const
{
    return std::make_unique<InterpolationPrismLinear>(*this);
}

Eigen::VectorXd InterpolationPrismLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd InterpolationPrismLinear::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationPrismLinear::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationPrismLinear::GetNumNodes() const
{
    return 6;
}

const Shape& InterpolationPrismLinear::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 3, 1> InterpolationPrismLinear::LocalCoords(int rNodeIndex)
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

Eigen::Matrix<double, 6, 1> InterpolationPrismLinear::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
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

Eigen::Matrix<double, 6, 3> InterpolationPrismLinear::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
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

std::vector<int> InterpolationPrismLinear::EdgeNodeIds(int edgeIndex) const
{
    switch (edgeIndex)
    {
    case 0:
        return {0, 1};
    case 1:
        return {1, 2};
    case 2:
        return {2, 0};
    case 3:
        return {3, 4};
    case 4:
        return {4, 5};
    case 5:
        return {5, 3};
    case 6:
        return {0, 3};
    case 7:
        return {1, 4};
    case 8:
        return {2, 5};
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "edge index out of range (0..8)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationPrismLinear::EdgeInterpolation(int /* edgeIndex*/) const
{
    return std::make_unique<InterpolationTrussLinear>();
}

std::vector<int> InterpolationPrismLinear::FaceNodeIds(int faceIndex) const
{
    switch (faceIndex)
    {
    case 0:
        return {0, 2, 1};
    case 1:
        return {3, 4, 5};
    case 2:
        return {0, 1, 4, 3};
    case 3:
        return {1, 2, 5, 4};
    case 4:
        return {0, 3, 5, 2};
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "edge index out of range (0..4)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationPrismLinear::FaceInterpolation(int faceIndex) const
{
    switch (faceIndex)
    {
    case 0:
    case 1:
        return std::make_unique<InterpolationTriangleLinear>();
    case 2:
    case 3:
    case 4:
        return std::make_unique<InterpolationQuadLinear>();
    default:
        throw Exception(__PRETTY_FUNCTION__, "Face index out of range (0-5).");
    }
}
