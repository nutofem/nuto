#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTrussQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTrussQuadratic>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 3, 1> ShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 3, 1> DerivativeShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTrussOrder2(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsTrussOrder2(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTrussOrder2(nodeId);
    }

    int GetNumNodes() const override
    {
        return 3;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Line mShape;
};
} /* NuTo */
