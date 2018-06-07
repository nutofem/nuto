#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTrussLinear : public InterpolationSimple
{
public:
    static Eigen::Matrix<double, 1, 1> LocalCoords(int rNodeIndex);

    static Eigen::Matrix<double, 2, 1> ShapeFunctions(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 2, 1> DerivativeShapeFunctions();

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    Eigen::MatrixXd GetDerivativeShapeFunctions(const NaturalCoords&) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

    std::vector<int> EdgeNodeIds(int) const override;

    int NumEdges() const override
    {
        return 1;
    }

    int NumFaces() const override
    {
        return 0;
    }

    virtual std::unique_ptr<InterpolationSimple> EdgeInterpolation(int /* edgeIndex*/) const override;

    virtual std::vector<int> FaceNodeIds(int /* faceIndex */) const override;

    virtual std::unique_ptr<InterpolationSimple> FaceInterpolation(int /* faceIndex*/) const override;

private:
    Line mShape;
};
} /* NuTo */
