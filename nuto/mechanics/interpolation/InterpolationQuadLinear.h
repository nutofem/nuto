#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Quadrilateral.h"

namespace NuTo
{
class InterpolationQuadLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationQuadLinear>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 4, 1> ShapeFunctionsQuadOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 4, 2> DerivativeShapeFunctionsQuadOrder1(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsQuadOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsQuadOrder1(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesQuadOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 4;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Quadrilateral mShape;
};
} /* NuTo */
