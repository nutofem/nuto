#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Triangle.h"

namespace NuTo
{
class InterpolationTriangleQuadratic : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder2(int rNodeIndex);

    static Eigen::Matrix<double, 6, 1> ShapeFunctionsTriangleOrder2(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 6, 2> DerivativeShapeFunctionsTriangleOrder2(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Triangle mShape;
};
} /* NuTo */
