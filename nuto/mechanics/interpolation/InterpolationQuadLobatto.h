#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/math/shapes/Quadrilateral.h"
#include "nuto/math/NaturalCoordinateMemoizer.h"

namespace NuTo
{
class InterpolationQuadLobatto : public InterpolationSimple
{
public:
    static Eigen::MatrixXd LocalCoords(int nodeId, const Eigen::VectorXd& nodes);

    static Eigen::VectorXd ShapeFunctions(const Eigen::Vector2d x, const Eigen::VectorXd& nodes);

    static Eigen::MatrixXd DerivativeShapeFunctions(const Eigen::Vector2d x, const Eigen::VectorXd& nodes);
    InterpolationQuadLobatto(int order);

    InterpolationQuadLobatto(const InterpolationQuadLobatto& other);

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    Eigen::MatrixXd GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Eigen::VectorXd mNodes;
    Quadrilateral mShape;
    NaturalCoordinateMemoizerMap<Eigen::VectorXd, NaturalCoords> mShapeFunctionMemo;
    NaturalCoordinateMemoizerMap<Eigen::MatrixXd, NaturalCoords> mShapeFunctionDerivativesMemo;
};
} /* NuTo */
