#include "nuto/base/Exception.h"
#include "InterpolationTetrahedronQuadratic.h"

using namespace NuTo;

std::unique_ptr<InterpolationSimple> InterpolationTetrahedronQuadratic::Clone() const
{
    return std::make_unique<InterpolationTetrahedronQuadratic>(*this);
}

Eigen::VectorXd InterpolationTetrahedronQuadratic::GetShapeFunctions(const NaturalCoords& naturalIpCoords) const
{
    return ShapeFunctions(naturalIpCoords);
}

Eigen::MatrixXd
InterpolationTetrahedronQuadratic::GetDerivativeShapeFunctions(const NaturalCoords& coords) const
{
    return DerivativeShapeFunctions(coords);
}

NaturalCoords InterpolationTetrahedronQuadratic::GetLocalCoords(int nodeId) const
{
    return LocalCoords(nodeId);
}

int InterpolationTetrahedronQuadratic::GetNumNodes() const
{
    return 10;
}

const Shape& InterpolationTetrahedronQuadratic::GetShape() const
{
    return mShape;
}

Eigen::Matrix<double, 3, 1> InterpolationTetrahedronQuadratic::LocalCoords(int rNodeIndex)
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
    case 4:
        return Eigen::Vector3d(0.5, 0.0, 0.0);
    case 5:
        return Eigen::Vector3d(0.5, 0.5, 0.0);
    case 6:
        return Eigen::Vector3d(0.0, 0.5, 0.0);
    case 7:
        return Eigen::Vector3d(0.0, 0.0, 0.5);
    case 8:
        return Eigen::Vector3d(0.0, 0.5, 0.5);
    case 9:
        return Eigen::Vector3d(0.5, 0.0, 0.5);
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..9)");
    }
}

Eigen::Matrix<double, 10, 1> InterpolationTetrahedronQuadratic::ShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 10, 1> shapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];
    double t = rCoordinates[2];

    // node 0 (0,0,0)
    shapeFunctions[0] = 1. - 3. * (r + s + t) + 2. * (r * r + s * s + t * t) + 4. * (r * s + r * t + s * t);

    // node 1 (1,0,0)
    shapeFunctions[1] = -r + 2. * r * r;

    // node 2 (0,1,0)
    shapeFunctions[2] = -s + 2. * s * s;

    // node 3 (0,0,1)
    shapeFunctions[3] = -t + 2. * t * t;

    // node 4 (0.5,0,0)
    shapeFunctions[4] = 4. * r * (1. - r - s - t);

    // node 5 (0.5,0.5,0)
    shapeFunctions[5] = 4. * r * s;

    // node 6 (0,0.5,0)
    shapeFunctions[6] = 4. * s * (1. - r - s - t);

    // node 7 (0,0,0.5)
    shapeFunctions[7] = 4. * t * (1. - r - s - t);

    // node 8 (0,0.5,0.5)
    shapeFunctions[8] = 4. * s * t;

    // node 9 (0.5,0,0.5)
    shapeFunctions[9] = 4. * r * t;

    return shapeFunctions;
}

Eigen::Matrix<double, 10, 3>
InterpolationTetrahedronQuadratic::DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 10, 3> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];
    double t = rCoordinates[2];

    derivativeShapeFunctions(0, 0) = -3. + 4. * (r + s + t);
    derivativeShapeFunctions(0, 1) = -3. + 4. * (r + s + t);
    derivativeShapeFunctions(0, 2) = -3. + 4. * (r + s + t);

    derivativeShapeFunctions(1, 0) = -1. + 4. * r;
    derivativeShapeFunctions(1, 1) = 0;
    derivativeShapeFunctions(1, 2) = 0;

    derivativeShapeFunctions(2, 0) = 0;
    derivativeShapeFunctions(2, 1) = -1. + 4. * s;
    derivativeShapeFunctions(2, 2) = 0;

    derivativeShapeFunctions(3, 0) = 0;
    derivativeShapeFunctions(3, 1) = 0;
    derivativeShapeFunctions(3, 2) = -1. + 4. * t;

    derivativeShapeFunctions(4, 0) = 4. - 8. * r - 4. * s - 4. * t;
    derivativeShapeFunctions(4, 1) = -4. * r;
    derivativeShapeFunctions(4, 2) = -4. * r;

    derivativeShapeFunctions(5, 0) = 4. * s;
    derivativeShapeFunctions(5, 1) = 4. * r;
    derivativeShapeFunctions(5, 2) = 0.;

    derivativeShapeFunctions(6, 0) = -4. * s;
    derivativeShapeFunctions(6, 1) = 4. - 8. * s - 4. * r - 4. * t;
    derivativeShapeFunctions(6, 2) = -4. * s;

    derivativeShapeFunctions(7, 0) = -4. * t;
    derivativeShapeFunctions(7, 1) = -4. * t;
    derivativeShapeFunctions(7, 2) = 4. - 8. * t - 4. * r - 4. * s;

    derivativeShapeFunctions(8, 0) = 0.;
    derivativeShapeFunctions(8, 1) = 4. * t;
    derivativeShapeFunctions(8, 2) = 4. * s;

    derivativeShapeFunctions(9, 0) = 4. * t;
    derivativeShapeFunctions(9, 1) = 0;
    derivativeShapeFunctions(9, 2) = 4. * r;

    return derivativeShapeFunctions;
}
