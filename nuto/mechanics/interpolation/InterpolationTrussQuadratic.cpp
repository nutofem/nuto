#include "nuto/base/Exception.h"
#include "InterpolationTrussQuadratic.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTrussQuadratic::Clone() const
{
    return std::make_unique<InterpolationTrussQuadratic>(*this);
}

Eigen::VectorXd InterpolationTrussQuadratic::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd InterpolationTrussQuadratic::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationTrussQuadratic::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationTrussQuadratic::GetNumNodes() const
{
    return 3;
}

const Shape& InterpolationTrussQuadratic::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 1, 1> InterpolationTrussQuadratic::LocalCoords(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1:
        return Eigen::Matrix<double, 1, 1>::Constant(0.);
    case 2:
        return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..2)");
    }
}

Eigen::Matrix<double, 3, 1> InterpolationTrussQuadratic::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 1> shapeFunctions;
    shapeFunctions[0] = 0.5 * (1. - rCoordinates[0]) - 0.5 * (1. - rCoordinates[0] * rCoordinates[0]);
    shapeFunctions[1] = 1. - rCoordinates[0] * rCoordinates[0];
    shapeFunctions[2] = 0.5 * (1. + rCoordinates[0]) - 0.5 * (1. - rCoordinates[0] * rCoordinates[0]);
    return shapeFunctions;
}

Eigen::Matrix<double, 3, 1> InterpolationTrussQuadratic::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 1> derivativeShapeFunctions;
    derivativeShapeFunctions[0] = -0.5 + rCoordinates[0];
    derivativeShapeFunctions[1] = -2.0 * rCoordinates[0];
    derivativeShapeFunctions[2] = 0.5 + rCoordinates[0];
    return derivativeShapeFunctions;
}

std::vector<int> InterpolationTrussQuadratic::EdgeNodeIds(int edgeIndex) const
{
    switch (edgeIndex)
    {
    case 0:
        return {0, 1, 2};
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "edge index out of range (0)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationTrussQuadratic::EdgeInterpolation(int /* edgeIndex*/) const
{
    return std::make_unique<InterpolationTrussQuadratic>();
}

std::vector<int> InterpolationTrussQuadratic::FaceNodeIds(int /* faceIndex */) const
{
    throw Exception(__PRETTY_FUNCTION__, "Truss interpolation has no faces");
}

std::unique_ptr<InterpolationSimple> InterpolationTrussQuadratic::FaceInterpolation(int /* faceIndex*/) const
{
    throw Exception(__PRETTY_FUNCTION__, "Truss interpolation has no faces");
}
