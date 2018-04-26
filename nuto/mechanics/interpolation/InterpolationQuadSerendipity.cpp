#include "nuto/base/Exception.h"
#include "InterpolationQuadSerendipity.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationQuadQuadratic::Clone() const
{
    return std::make_unique<InterpolationQuadQuadratic>(*this);
}

Eigen::VectorXd InterpolationQuadQuadratic::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

DerivativeShapeFunctionsNatural
InterpolationQuadQuadratic::GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return DerivativeShapeFunctions(naturalIpCoords);
}

NaturalCoords InterpolationQuadQuadratic::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationQuadQuadratic::GetNumNodes() const
{
    return 8;
}

const Shape& InterpolationQuadQuadratic::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 2, 1> InterpolationQuadQuadratic::LocalCoords(int rNodeIndex)
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
    case 4:
        return Eigen::Vector2d(0.0, -1.0);
    case 5:
        return Eigen::Vector2d(1.0, 0.0);
    case 6:
        return Eigen::Vector2d(0.0, 1.0);
    case 7:
        return Eigen::Vector2d(-1.0, 0.0);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..7)");
    }
}

Eigen::Matrix<double, 8, 1> InterpolationQuadQuadratic::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 8, 1> shapeFunctions;
    double r = rCoordinates(0);
    double s = rCoordinates(1);

    shapeFunctions[0] = -0.25 * (r - 1) * (s - 1) * (r + s + 1);
    shapeFunctions[1] = 0.25 * (r + 1) * (s - 1) * (-r + s + 1);
    shapeFunctions[2] = 0.25 * (r + 1) * (s + 1) * (r + s - 1);
    shapeFunctions[3] = 0.25 * (r - 1) * (s + 1) * (r - s + 1);
    shapeFunctions[4] = 0.5 * (r * r - 1) * (s - 1);
    shapeFunctions[5] = -0.5 * (r + 1) * (s * s - 1);
    shapeFunctions[6] = -0.5 * (r * r - 1) * (s + 1);
    shapeFunctions[7] = 0.5 * (r - 1) * (s * s - 1);

    return shapeFunctions;
}

Eigen::Matrix<double, 8, 2> InterpolationQuadQuadratic::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 8, 2> derivativeShapeFunctions;
    double r = rCoordinates(0);
    double s = rCoordinates(1);

    derivativeShapeFunctions(0, 0) = (2 * r + s) * (-0.25 * s + 0.25);
    derivativeShapeFunctions(0, 1) = (-0.25 * r + 0.25) * (r + 2 * s);

    derivativeShapeFunctions(1, 0) = (-0.5 * r + 0.25 * s) * (s - 1);
    derivativeShapeFunctions(1, 1) = 0.25 * (-r + 2 * s) * (r + 1);

    derivativeShapeFunctions(2, 0) = 0.25 * (2 * r + s) * (s + 1);
    derivativeShapeFunctions(2, 1) = 0.25 * (r + 1) * (r + 2 * s);

    derivativeShapeFunctions(3, 0) = 0.25 * (2 * r - s) * (s + 1);
    derivativeShapeFunctions(3, 1) = (0.25 * r - 0.5 * s) * (r - 1);

    derivativeShapeFunctions(4, 0) = r * (s - 1);
    derivativeShapeFunctions(4, 1) = 0.5 * r * r - 0.5;

    derivativeShapeFunctions(5, 0) = -0.5 * s * s + 0.5;
    derivativeShapeFunctions(5, 1) = -s * (r + 1);

    derivativeShapeFunctions(6, 0) = -r * (s + 1);
    derivativeShapeFunctions(6, 1) = -0.5 * r * r + 0.5;

    derivativeShapeFunctions(7, 0) = 0.5 * s * s - 0.5;
    derivativeShapeFunctions(7, 1) = s * (r - 1);

    return derivativeShapeFunctions;
}
