#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Prism.h"

namespace NuTo
{
class InterpolationPrismQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesPrismOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 18, 1> ShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 18, 3> DerivativeShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Prism mShape;
};
} /* NuTo */
