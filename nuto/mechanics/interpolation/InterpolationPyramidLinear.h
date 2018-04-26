#pragma once
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Pyramid.h"

namespace NuTo
{
class InterpolationPyramidLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesPyramidOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 5, 1> ShapeFunctionsPyramidOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 5, 3> DerivativeShapeFunctionsPyramidOrder1(const Eigen::VectorXd& rCoordinates);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Pyramid mShape;
};
} /* NuTo */
