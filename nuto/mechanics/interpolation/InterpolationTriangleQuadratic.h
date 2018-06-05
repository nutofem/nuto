#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Triangle.h"
#include <vector>

namespace NuTo
{
class InterpolationTriangleQuadratic : public InterpolationSimple
{
public:
    static Eigen::Matrix<double, 2, 1> LocalCoords(int rNodeIndex);

    static Eigen::Matrix<double, 6, 1> ShapeFunctions(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 6, 2> DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates);

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    Eigen::MatrixXd GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

    std::vector<int> EdgeNodeIds(int edgeIndex) const override;

    std::unique_ptr<InterpolationSimple> EdgeInterpolation(int edgeIndex) const override;

    int NumEdges() const override
    {
        return 3;
    }

private:
    Triangle mShape;
};
} /* NuTo */
