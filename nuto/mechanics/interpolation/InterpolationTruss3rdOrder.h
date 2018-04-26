#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Line.h"

namespace NuTo
{
class InterpolationTruss3rdOrder : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder3(int rNodeIndex);

    static Eigen::Matrix<double, 4, 1> ShapeFunctionsTrussOrder3(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 4, 1> DerivativeShapeFunctionsTrussOrder3(const Eigen::VectorXd& rCoordinates);
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
