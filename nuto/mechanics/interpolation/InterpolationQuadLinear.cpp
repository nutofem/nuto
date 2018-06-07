#include "nuto/base/Exception.h"
#include "InterpolationQuadLinear.h"
#include "InterpolationTrussLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationQuadLinear::Clone() const
{
    return std::make_unique<InterpolationQuadLinear>(*this);
}

Eigen::VectorXd InterpolationQuadLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd InterpolationQuadLinear::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationQuadLinear::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationQuadLinear::GetNumNodes() const
{
    return 4;
}

const Shape& InterpolationQuadLinear::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 2, 1> InterpolationQuadLinear::LocalCoords(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector2d(-1.0, -1.0);
    case 1:
        return Eigen::Vector2d(1.0, -1.0);
    case 2:
        return Eigen::Vector2d(1.0, 1.0);
    case 3:
        return Eigen::Vector2d(-1.0, 1.0);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..3)");
    }
}

Eigen::Matrix<double, 4, 1> InterpolationQuadLinear::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> shapeFunctions;
    shapeFunctions[0] = 0.25 * (1. - rCoordinates(0)) * (1. - rCoordinates(1));
    shapeFunctions[1] = 0.25 * (1. + rCoordinates(0)) * (1. - rCoordinates(1));
    shapeFunctions[2] = 0.25 * (1. + rCoordinates(0)) * (1. + rCoordinates(1));
    shapeFunctions[3] = 0.25 * (1. - rCoordinates(0)) * (1. + rCoordinates(1));
    return shapeFunctions;
}

Eigen::Matrix<double, 4, 2> InterpolationQuadLinear::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 2> derivativeShapeFunctions;
    derivativeShapeFunctions(0, 0) = -0.25 * (1. - rCoordinates(1));
    derivativeShapeFunctions(0, 1) = -0.25 * (1. - rCoordinates(0));

    derivativeShapeFunctions(1, 0) = +0.25 * (1. - rCoordinates(1));
    derivativeShapeFunctions(1, 1) = -0.25 * (1. + rCoordinates(0));

    derivativeShapeFunctions(2, 0) = +0.25 * (1. + rCoordinates(1));
    derivativeShapeFunctions(2, 1) = +0.25 * (1. + rCoordinates(0));

    derivativeShapeFunctions(3, 0) = -0.25 * (1. + rCoordinates(1));
    derivativeShapeFunctions(3, 1) = +0.25 * (1. - rCoordinates(0));

    return derivativeShapeFunctions;
}

std::vector<int> InterpolationQuadLinear::EdgeNodeIds(int edgeIndex) const
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
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "edge index out of range (0..3)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationQuadLinear::EdgeInterpolation(int /* edgeIndex*/) const
{
    return std::make_unique<InterpolationTrussLinear>();
}

std::vector<int> InterpolationQuadLinear::FaceNodeIds(int faceIndex) const
{
    switch (faceIndex)
    {
    case 0:
        return {0, 1, 2, 3};
    default:
        throw Exception(__PRETTY_FUNCTION__, "Face index out of range (0).");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationQuadLinear::FaceInterpolation(int /* faceIndex*/) const
{
    return std::make_unique<InterpolationQuadLinear>();
}
