#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Tetrahedron.h"

namespace NuTo
{
class InterpolationTetrahedronLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesTetrahedronOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 4, 1> ShapeFunctionsTetrahedronOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 4, 3> DerivativeShapeFunctionsTetrahedronOrder1();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override;

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords&) const override;

    NaturalCoords GetLocalCoords(int nodeId) const override;

    int GetNumNodes() const override;

    const Shape& GetShape() const override;

private:
    Tetrahedron mShape;
};
} /* NuTo */
