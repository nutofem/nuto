#pragma once

#include "mechanics/cell/Jacobian.h"
#include "mechanics/interpolation/ElementInterpolationBase.h"
#include "mechanics/nodes/DofVector.h"
#include "mechanics/nodes/DofVector.h"

namespace NuTo
{

template <int TDim>
class PDE_Element
{
public:
    PDE_Element(const ElementInterpolationBase& coordElement, DofContainer<ElementInterpolationBase*> elements)
        : mCoordinateElement(coordElement)
        , mElements(elements)
    {
    }

    NodeValues ExtractNodeValues(const DofType& dofType) const
    {
        return mElements[dofType]->ExtractNodeValues();
    }

    DofVector<int> DofNumbering()
    {
        return DofVector<int>();
    }

    Jacobian<TDim> ComputeJacobian(const NaturalCoords& naturalIPCoords)const
    {
        return Jacobian<TDim>(mCoordinateElement.ExtractNodeValues(),
                        mCoordinateElement.GetDerivativeShapeFunctions(naturalIPCoords));
    }

    const ElementInterpolationBase* GetInterpolation(const DofType& dofType)const
    {
        return mElements[dofType];
    }


private:
    const ElementInterpolationBase& mCoordinateElement;
    DofContainer<ElementInterpolationBase*> mElements;
};
} /* NuTo */
