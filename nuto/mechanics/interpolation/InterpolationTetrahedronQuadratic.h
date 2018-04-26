#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Tetrahedron.h"

namespace NuTo
{
class InterpolationTetrahedronQuadratic : public InterpolationSimple
{
public:
    static Eigen::Matrix<double, 3, 1> LocalCoords(int rNodeIndex);

    static Eigen::Matrix<double, 10, 1> ShapeFunctions(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 10, 3> DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates);

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& coords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Tetrahedron mShape;
};
} /* NuTo */
