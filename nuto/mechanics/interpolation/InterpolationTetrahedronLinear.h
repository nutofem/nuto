#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/InterpolationSimple.h"
#include "nuto/math/shapes/Tetrahedron.h"

namespace NuTo
{
class InterpolationTetrahedronLinear : public InterpolationSimple
{
public:
    std::unique_ptr<InterpolationSimple> Clone() const override
    {
        return std::make_unique<InterpolationTetrahedronLinear>(*this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static Eigen::Matrix<double, 3, 1> NodeCoordinatesTetrahedronOrder1(int rNodeIndex);

    static Eigen::Matrix<double, 4, 1> ShapeFunctionsTetrahedronOrder1(const Eigen::VectorXd& rCoordinates);

    static Eigen::Matrix<double, 4, 3> DerivativeShapeFunctionsTetrahedronOrder1();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ShapeFunctions GetShapeFunctions(const NaturalCoords& naturalIpCoords) const override
    {
        return ShapeFunctionsTetrahedronOrder1(naturalIpCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords&) const override
    {
        return DerivativeShapeFunctionsTetrahedronOrder1();
    }

    NaturalCoords GetLocalCoords(int nodeId) const override
    {
        return NodeCoordinatesTetrahedronOrder1(nodeId);
    }

    int GetNumNodes() const override
    {
        return 4;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    Tetrahedron mShape;
};
} /* NuTo */
