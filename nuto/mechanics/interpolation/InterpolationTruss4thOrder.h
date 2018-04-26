#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTruss4thOrder : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder4(int rNodeIndex);

    static Eigen::Matrix<double, 5, 1> ShapeFunctionsTrussOrder4(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 5, 1> DerivativeShapeFunctionsTrussOrder4(const Eigen::VectorXd& rCoordinates);

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
