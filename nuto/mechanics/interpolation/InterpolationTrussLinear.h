#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTrussLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTrussLinear>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 2, 1> ShapeFunctionsTrussOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 2, 1> DerivativeShapeFunctionsTrussOrder1();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTrussOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords&) const override
    {
        return DerivativeShapeFunctionsTrussOrder1();
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTrussOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 2;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Line mShape;
};
} /* NuTo */
