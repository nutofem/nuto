#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTruss3rdOrder : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTruss3rdOrder>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder3(int rNodeIndex);

    static Eigen::Matrix<double, 4, 1> ShapeFunctionsTrussOrder3(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 4, 1> DerivativeShapeFunctionsTrussOrder3(const Eigen::VectorXd& rCoordinates);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTrussOrder3(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return DerivativeShapeFunctionsTrussOrder3(naturalIpCoords);
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTrussOrder3(nodeId);
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
    Line mShape;
};
} /* NuTo */
