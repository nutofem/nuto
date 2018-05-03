#include "nuto/base/Exception.h"
#include "InterpolationTetrahedronLinear.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTetrahedronLinear::Clone() const
{
    return std::make_unique<InterpolationTetrahedronLinear>(*this);
}

Eigen::VectorXd InterpolationTetrahedronLinear::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd InterpolationTetrahedronLinear::GetDerivativeShapeFunctions(const NaturalCoords&) const
{
    return DerivativeShapeFunctions();
}

NaturalCoords InterpolationTetrahedronLinear::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationTetrahedronLinear::GetNumNodes() const
{
    return 4;
}

const Shape& InterpolationTetrahedronLinear::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 3, 1> InterpolationTetrahedronLinear::LocalCoords(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector3d(0.0, 0.0, 0.0);
    case 1:
        return Eigen::Vector3d(1.0, 0.0, 0.0);
    case 2:
        return Eigen::Vector3d(0.0, 1.0, 0.0);
    case 3:
        return Eigen::Vector3d(0.0, 0.0, 1.0);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..4)");
    }
}

Eigen::Matrix<double, 4, 1> InterpolationTetrahedronLinear::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> shapeFunctions;
    shapeFunctions[0] = 1. - rCoordinates.sum();
    shapeFunctions[1] = rCoordinates(0);
    shapeFunctions[2] = rCoordinates(1);
    shapeFunctions[3] = rCoordinates(2);
    return shapeFunctions;
}

Eigen::Matrix<double, 4, 3> InterpolationTetrahedronLinear::DerivativeShapeFunctions()
{
    Eigen::Matrix<double, 4, 3> derivativeShapeFunctions;
    derivativeShapeFunctions(0, 0) = -1;
    derivativeShapeFunctions(0, 1) = -1;
    derivativeShapeFunctions(0, 2) = -1;

    derivativeShapeFunctions(1, 0) = 1;
    derivativeShapeFunctions(1, 1) = 0;
    derivativeShapeFunctions(1, 2) = 0;

    derivativeShapeFunctions(2, 0) = 0;
    derivativeShapeFunctions(2, 1) = 1;
    derivativeShapeFunctions(2, 2) = 0;

    derivativeShapeFunctions(3, 0) = 0;
    derivativeShapeFunctions(3, 1) = 0;
    derivativeShapeFunctions(3, 2) = 1;

    return derivativeShapeFunctions;
}
