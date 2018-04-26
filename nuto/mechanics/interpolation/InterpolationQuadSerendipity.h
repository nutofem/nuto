#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Quadrilateral.h"

namespace NuTo
{
class InterpolationQuadQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationQuadQuadratic>(*this);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 8, 1> ShapeFunctionsQuadOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 8, 2> DerivativeShapeFunctionsQuadOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsQuadOrder2(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsQuadOrder2(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesQuadOrder2(nodeId);
    }

    int GetNumNodes() const override
    {
        return 8;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Quadrilateral mShape;
};
} /* NuTo */
