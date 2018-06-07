#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/math/shapes/Quadrilateral.h"

namespace NuTo
{
class InterpolationQuadLobatto : public InterpolationSimple
{
public:
    static Eigen::MatrixXd LocalCoords(int nodeId, const Eigen::VectorXd& nodes);

    static Eigen::VectorXd ShapeFunctions(const Eigen::Vector2d x, const Eigen::VectorXd& nodes);

    static Eigen::MatrixXd DerivativeShapeFunctions(const Eigen::Vector2d x, const Eigen::VectorXd& nodes);
    InterpolationQuadLobatto(int order);

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    Eigen::MatrixXd GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

    std::vector<int> EdgeNodeIds(int edgeIndex) const override;

    int NumEdges() const override
    {
        return 4;
    }

    int NumFaces() const override
    {
        return 1;
    }

    virtual std::unique_ptr<InterpolationSimple> EdgeInterpolation(int /* edgeIndex*/) const override;

    virtual std::vector<int> FaceNodeIds(int /* faceIndex */) const override;

    virtual std::unique_ptr<InterpolationSimple> FaceInterpolation(int /* faceIndex*/) const override;

private:
    Eigen::VectorXd mNodes;
    Quadrilateral mShape;
};
} /* NuTo */
