#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTrussLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 2, 1> ShapeFunctionsTrussOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 2, 1> DerivativeShapeFunctionsTrussOrder1();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords&) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Line mShape;
};
} /* NuTo */
