#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Pyramid.h"

namespace NuTo
{
class InterpolationPyramidLinear : public InterpolationSimple
{
public:
    static Eigen::Matrix<double, 3, 1> LocalCoords(int rNodeIndex);

    static Eigen::Matrix<double, 5, 1> ShapeFunctions(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 5, 3> DerivativeShapeFunctions(const Eigen::VectorXd& rCoordinates);

    std::unique_ptr<InterpolationSimple> Clone() const override;

    Eigen::VectorXd GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Pyramid mShape;
};
} /* NuTo */
