#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/Legendre.h"
#include "nuto/math/shapes/Line.h"
#include "nuto/base/Exception.h"

namespace NuTo
{
class InterpolationTrussLobatto : public InterpolationSimple
{
public:

    static Eigen::VectorXd LocalCoords(int order);

    static Eigen::VectorXd BarycentricWeights(const Eigen::VectorXd& nodes);

    static Eigen::VectorXd ShapeFunctions(const double x, const Eigen::VectorXd& nodes);

    static Eigen::VectorXd DerivativeShapeFunctions(const double x, const Eigen::VectorXd& nodes);

    InterpolationTrussLobatto(int order);

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    Eigen::MatrixXd GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Eigen::VectorXd mNodes;
    Line mShape;
};
} /* NuTo */
