#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTruss4thOrder : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTruss4thOrder>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder4(int rNodeIndex);

    static Eigen::Matrix<double, 5, 1> ShapeFunctionsTrussOrder4(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 5, 1> DerivativeShapeFunctionsTrussOrder4(const Eigen::VectorXd& rCoordinates);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTrussOrder4(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsTrussOrder4(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTrussOrder4(nodeId);
    }

    int GetNumNodes() const override
    {
        return 5;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Line mShape;
};
} /* NuTo */
