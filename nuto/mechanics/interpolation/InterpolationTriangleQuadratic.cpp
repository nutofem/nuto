#include "nuto/base/Exception.h"
#include "InterpolationTriangleQuadratic.h"
#include "InterpolationTrussQuadratic.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTriangleQuadratic::Clone() const
{
    return std::make_unique<InterpolationTriangleQuadratic>(*this);
}

Eigen::VectorXd InterpolationTriangleQuadratic::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd InterpolationTriangleQuadratic::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationTriangleQuadratic::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationTriangleQuadratic::GetNumNodes() const
{
    return 6;
}

const Shape& InterpolationTriangleQuadratic::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 2, 1> InterpolationTriangleQuadratic::LocalCoords(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector2d(0.0, 0.0);
    case 1:
        return Eigen::Vector2d(1.0, 0.0);
    case 2:
        return Eigen::Vector2d(0.0, 1.0);
    case 3:
        return Eigen::Vector2d(0.5, 0.0);
    case 4:
        return Eigen::Vector2d(0.5, 0.5);
    case 5:
        return Eigen::Vector2d(0.0, 0.5);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..5)");
    }
}

Eigen::Matrix<double, 6, 1> InterpolationTriangleQuadratic::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 6, 1> shapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];

    shapeFunctions[0] = 2. * (r * r + s * s) + 4. * r * s - 3. * (r + s) + 1.;
    shapeFunctions[1] = 2. * r * r - r;
    shapeFunctions[2] = 2. * s * s - s;
    shapeFunctions[3] = -4. * r * (r + s - 1.);
    shapeFunctions[4] = 4. * r * s;
    shapeFunctions[5] = -4. * s * (s + r - 1.);
    return shapeFunctions;
}

Eigen::Matrix<double, 6, 2>
InterpolationTriangleQuadratic::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 6, 2> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];

    derivativeShapeFunctions(0, 0) = 4. * (r + s) - 3.;
    derivativeShapeFunctions(0, 1) = 4. * (r + s) - 3.;

    derivativeShapeFunctions(1, 0) = 4. * r - 1.;
    derivativeShapeFunctions(1, 1) = 0.;

    derivativeShapeFunctions(2, 0) = 0.;
    derivativeShapeFunctions(2, 1) = 4. * s - 1.;

    derivativeShapeFunctions(3, 0) = -8. * r - 4. * s + 4.;
    derivativeShapeFunctions(3, 1) = -4. * r;

    derivativeShapeFunctions(4, 0) = 4. * s;
    derivativeShapeFunctions(4, 1) = 4. * r;

    derivativeShapeFunctions(5, 0) = -4. * s;
    derivativeShapeFunctions(5, 1) = -8. * s - 4. * r + 4.;
    return derivativeShapeFunctions;
}

std::vector<int> InterpolationTriangleQuadratic::EdgeNodeIds(int edgeIndex) const
{
    switch (edgeIndex)
    {
    case 0:
        return {0, 3, 1};
    case 1:
        return {1, 4, 2};
    case 2:
        return {2, 5, 0};
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "edge index out of range (0..2)");
    }
}

std::unique_ptr<InterpolationSimple> InterpolationTriangleQuadratic::EdgeInterpolation(int /* edgeIndex */) const
{
    return std::make_unique<InterpolationTrussQuadratic>();
}
