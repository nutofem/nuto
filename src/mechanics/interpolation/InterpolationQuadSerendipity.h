#pragma once
#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/elements/ElementShapeFunctions.h"

namespace NuTo
{
class InterpolationQuadSerendipity : public Interpolation
{
public:
    InterpolationQuadSerendipity(int rDofDimension)
        : mDofDimension(rDofDimension)
    {
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        return ShapeFunctions2D::ShapeFunctionsQuadOrder2(rNaturalIPCoords);
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder2(rNaturalIPCoords);
    }

    NaturalCoords GetLocalCoords(int rNodeId) const override
    {
        return ShapeFunctions2D::NodeCoordinatesQuadOrder2(rNodeId);
    }

    int GetNumNodes() const override
    {
        return 8;
    }

    int GetDofDimension() const override
    {
        return mDofDimension;
    }

private:
    int mDofDimension;
};
} /* NuTo */
