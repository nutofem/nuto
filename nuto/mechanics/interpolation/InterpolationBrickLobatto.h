#pragma once

#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/math/shapes/Hexahedron.h"

namespace NuTo
{
class InterpolationBrickLobatto : public InterpolationSimple
{
public:

    static Eigen::MatrixXd LocalCoords(int nodeId, const Eigen::VectorXd& nodes);

    static Eigen::VectorXd ShapeFunctions(const Eigen::Vector3d x, const Eigen::VectorXd& nodes);

    static Eigen::MatrixXd DerivativeShapeFunctions(const Eigen::Vector3d x,
                                                                            const Eigen::VectorXd& nodes);

    InterpolationBrickLobatto(int order);

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Eigen::VectorXd mNodes;
    Hexahedron mShape;
};
} /* NuTo */
