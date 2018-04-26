#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTrussQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 3, 1> ShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 3, 1> DerivativeShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Line mShape;
};
} /* NuTo */
