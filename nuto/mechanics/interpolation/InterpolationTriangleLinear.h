#pragma once
#include <Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Triangle.h"

namespace NuTo
{
class InterpolationTriangleLinear : public InterpolationSimple
{
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 3, 1> ShapeFunctionsTriangleOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 3, 2> DerivativeShapeFunctionsTriangleOrder1(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<InterpolationSimple> Clone() const override;

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Triangle mShape;
};
} /* NuTo */
