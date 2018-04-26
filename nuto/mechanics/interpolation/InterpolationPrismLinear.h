#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Prism.h"

namespace NuTo
{
class InterpolationPrismLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationPrismLinear>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesPrismOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 6, 1> ShapeFunctionsPrismOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 6, 3> DerivativeShapeFunctionsPrismOrder1(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsPrismOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsPrismOrder1(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesPrismOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 6;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Prism mShape;
};
} /* NuTo */
