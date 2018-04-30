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

    static Eigen::VectorXd NodeCoordinatesTrussLobatto(int order);

    static Eigen::VectorXd BarycentricWeights(const Eigen::VectorXd& nodes);

    static Eigen::VectorXd ShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes);

    static Eigen::VectorXd DerivativeShapeFunctionsTrussLagrange(const double x, const Eigen::VectorXd& nodes);

    InterpolationTrussLobatto(int order);

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Eigen::VectorXd mNodes;
    Line mShape;
};
} /* NuTo */
