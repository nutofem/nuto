#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Prism.h"

namespace NuTo
{
class InterpolationPrismQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationPrismQuadratic>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesPrismOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 18, 1> ShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 18, 3> DerivativeShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsPrismOrder2(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsPrismOrder2(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesPrismOrder2(nodeId);
    }

    int GetNumNodes() const override
    {
        return 18;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Prism mShape;
};
} /* NuTo */
