#pragma once
#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/elements/ElementShapeFunctions.h"

namespace NuTo
{
class InterpolationQuadLinear : public Interpolation
{
public:
    InterpolationQuadLinear(int rDofDimension)
        : mDofDimension(rDofDimension)
    {
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        return ShapeFunctions2D::ShapeFunctionsQuadOrder1(rNaturalIPCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder1(rNaturalIPCoords);
    }

    NaturalCoords GetLocalCoords(int rNodeId) const override
    {
        return ShapeFunctions2D::NodeCoordinatesQuadOrder1(rNodeId);
    }

    int GetNumNodes() const override
    {
        return 4;
    }

    int GetDofDimension() const override
    {
        return mDofDimension;
    }

private:
    int mDofDimension;
};
} /* NuTo */
