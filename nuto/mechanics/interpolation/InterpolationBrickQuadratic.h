#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Hexahedron.h"

namespace NuTo
{
class InterpolationBrickQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 20, 1> ShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 20, 3> DerivativeShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Hexahedron mShape;
};
} /* NuTo */
