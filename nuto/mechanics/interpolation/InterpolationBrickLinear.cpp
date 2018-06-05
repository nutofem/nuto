#include "nuto/base/Exception.h"
#include "InterpolationBrickLinear.h"
#include "InterpolationQuadLinear.h"
#include "InterpolationTrussLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationBrickLinear::Clone() const
{
    return std::make_unique<InterpolationBrickLinear>(*this);
}

Eigen::VectorXd InterpolationBrickLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd InterpolationBrickLinear::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationBrickLinear::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationBrickLinear::GetNumNodes() const
{
    return 8;
}

const Shape& InterpolationBrickLinear::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 3, 1> InterpolationBrickLinear::LocalCoords(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0:
        return Eigen::Vector3d(-1., -1., -1.);
    case 1:
        return Eigen::Vector3d(1., -1., -1.);
    case 2:
        return Eigen::Vector3d(1., 1., -1.);
    case 3:
        return Eigen::Vector3d(-1., 1., -1.);
    case 4:
        return Eigen::Vector3d(-1., -1., 1.);
    case 5:
        return Eigen::Vector3d(1., -1., 1.);
    case 6:
        return Eigen::Vector3d(1., 1., 1.);
    case 7:
        return Eigen::Vector3d(-1., 1., 1.);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..7)");
    }
}

Eigen::Matrix<double, 8, 1> InterpolationBrickLinear::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 8, 1> shapeFunctions;

    double plus_r = 1.0 + rCoordinates[0];
    double plus_s = 1.0 + rCoordinates[1];
    double plus_t = 1.0 + rCoordinates[2];

    double minus_r = 1.0 - rCoordinates[0];
    double minus_s = 1.0 - rCoordinates[1];
    double minus_t = 1.0 - rCoordinates[2];

    shapeFunctions[0] = 0.125 * minus_r * minus_s * minus_t;
    shapeFunctions[1] = 0.125 * plus_r * minus_s * minus_t;
    shapeFunctions[2] = 0.125 * plus_r * plus_s * minus_t;
    shapeFunctions[3] = 0.125 * minus_r * plus_s * minus_t;
    shapeFunctions[4] = 0.125 * minus_r * minus_s * plus_t;
    shapeFunctions[5] = 0.125 * plus_r * minus_s * plus_t;
    shapeFunctions[6] = 0.125 * plus_r * plus_s * plus_t;
    shapeFunctions[7] = 0.125 * minus_r * plus_s * plus_t;

    return shapeFunctions;
}

Eigen::Matrix<double, 8, 3> InterpolationBrickLinear::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    double plus_r = 1.0 + rCoordinates[0];
    double plus_s = 1.0 + rCoordinates[1];
    double plus_t = 1.0 + rCoordinates[2];

    double minus_r = 1.0 - rCoordinates[0];
    double minus_s = 1.0 - rCoordinates[1];
    double minus_t = 1.0 - rCoordinates[2];

    Eigen::Matrix<double, 8, 3> derivativeShapeFunctions;

    derivativeShapeFunctions(0, 0) = -0.125 * minus_s * minus_t;
    derivativeShapeFunctions(0, 1) = -0.125 * minus_r * minus_t;
    derivativeShapeFunctions(0, 2) = -0.125 * minus_r * minus_s;

    derivativeShapeFunctions(1, 0) = 0.125 * minus_s * minus_t;
    derivativeShapeFunctions(1, 1) = -0.125 * plus_r * minus_t;
    derivativeShapeFunctions(1, 2) = -0.125 * plus_r * minus_s;

    derivativeShapeFunctions(2, 0) = 0.125 * plus_s * minus_t;
    derivativeShapeFunctions(2, 1) = 0.125 * plus_r * minus_t;
    derivativeShapeFunctions(2, 2) = -0.125 * plus_r * plus_s;

    derivativeShapeFunctions(3, 0) = -0.125 * plus_s * minus_t;
    derivativeShapeFunctions(3, 1) = 0.125 * minus_r * minus_t;
    derivativeShapeFunctions(3, 2) = -0.125 * minus_r * plus_s;

    derivativeShapeFunctions(4, 0) = -0.125 * minus_s * plus_t;
    derivativeShapeFunctions(4, 1) = -0.125 * minus_r * plus_t;
    derivativeShapeFunctions(4, 2) = 0.125 * minus_r * minus_s;

    derivativeShapeFunctions(5, 0) = 0.125 * minus_s * plus_t;
    derivativeShapeFunctions(5, 1) = -0.125 * plus_r * plus_t;
    derivativeShapeFunctions(5, 2) = 0.125 * plus_r * minus_s;

    derivativeShapeFunctions(6, 0) = 0.125 * plus_s * plus_t;
    derivativeShapeFunctions(6, 1) = 0.125 * plus_r * plus_t;
    derivativeShapeFunctions(6, 2) = 0.125 * plus_r * plus_s;

    derivativeShapeFunctions(7, 0) = -0.125 * plus_s * plus_t;
    derivativeShapeFunctions(7, 1) = 0.125 * minus_r * plus_t;
    derivativeShapeFunctions(7, 2) = 0.125 * minus_r * plus_s;

    return derivativeShapeFunctions;
}

std::vector<int> InterpolationBrickLinear::EdgeNodeIds(int edgeIndex) const
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
        return {4, 5};
    case 5:
        return {5, 6};
    case 6:
        return {6, 7};
    case 7:
        return {7, 4};
    case 8:
        return {0, 4};
    case 9:
        return {1, 5};
    case 10:
        return {2, 6};
    case 11:
        return {3, 7};
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "edge index out of range (0..11)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationBrickLinear::EdgeInterpolation(int /* faceIndex */) const
{
    return std::make_unique<InterpolationTrussLinear>();
}

std::vector<int> InterpolationBrickLinear::FaceNodeIds(int faceIndex) const
{
    switch (faceIndex)
    {
    case 0:
        return {0, 3, 2, 1};
    case 1:
        return {0, 1, 5, 4};
    case 2:
        return {4, 5, 6, 7};
    case 3:
        return {0, 4, 7, 3};
    case 4:
        return {2, 3, 7, 6};
    case 5:
        return {1, 2, 6, 5};
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "face index out of range (0..5)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationBrickLinear::FaceInterpolation(int /* faceIndex */) const
{
    return std::make_unique<InterpolationQuadLinear>();
}
