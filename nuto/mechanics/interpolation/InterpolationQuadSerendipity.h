#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Quadrilateral.h"

namespace NuTo
{
class InterpolationQuadQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 8, 1> ShapeFunctionsQuadOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 8, 2> DerivativeShapeFunctionsQuadOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Quadrilateral mShape;
};
} /* NuTo */
