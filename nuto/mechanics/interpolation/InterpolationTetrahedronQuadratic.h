#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Tetrahedron.h"

namespace NuTo
{
class InterpolationTetrahedronQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesTetrahedronOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 10, 1> ShapeFunctionsTetrahedronOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 10, 3> DerivativeShapeFunctionsTetrahedronOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& coords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Tetrahedron mShape;
};
} /* NuTo */
