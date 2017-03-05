#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/elements/ElementShapeFunctions.h"

namespace NuTo
{
class InterpolationTriangle : public Interpolation
{
public:
    InterpolationTriangle(eInterpolation rType, int rOrder, int rDofDimension)
        : mType(rType)
        , mOrder(rOrder)
        , mDofDimension(rDofDimension)
    {
    }

    ShapeFunctions GetShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        if (mType == eInterpolation::GAUSS && mOrder == 1)
            return ShapeFunctions2D::ShapeFunctionsTriangleOrder1(rNaturalIPCoords);
        throw;
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(const NaturalCoords& rNaturalIPCoords) const override
    {
        if (mType == eInterpolation::GAUSS && mOrder == 1)
            return ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder1(rNaturalIPCoords);
        throw;
    }

    NaturalCoords GetLocalCoords(int rNodeId) const override
    {
        if (mType == eInterpolation::GAUSS && mOrder == 1)
            return ShapeFunctions2D::NodeCoordinatesTriangleOrder1(rNodeId);
        throw;
    }

    int GetNumNodes() const override
    {
        if (mType == eInterpolation::GAUSS && mOrder == 1)
            return 3;
        throw;
    }

    int GetDofDimension() const override
    {
        return mDofDimension;
    }

private:
    eInterpolation mType;
    int mOrder;
    int mDofDimension;
};
} /* NuTo */
