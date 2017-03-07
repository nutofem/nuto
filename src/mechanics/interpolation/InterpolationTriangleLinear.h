#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/elements/ElementShapeFunctions.h"

namespace NuTo
{
class InterpolationTriangleLinear : public Interpolation
{
public:
    InterpolationTriangleLinear(int rDofDimension)
        : mDofDimension(rDofDimension)
    {
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        return ShapeFunctions2D::ShapeFunctionsTriangleOrder1(rNaturalIPCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        return ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder1(rNaturalIPCoords);
    }

    NaturalCoords GetLocalCoords(int rNodeId) const override
    {
        return ShapeFunctions2D::NodeCoordinatesTriangleOrder1(rNodeId);
    }

    int GetNumNodes() const override
    {
        return 3;
    }

    int GetDofDimension() const override
    {
        return mDofDimension;
    }

private:
    int mDofDimension;
};
} /* NuTo */
