#pragma once
#include "mechanics/elements/ElementInterface.h"
#include "mechanics/nodes/DofContainer.h"

namespace NuTo
{
class ElementCollection
{
public:
    ElementCollection(const ElementInterface& coordinateElement)
        : mCoordinateElement(coordinateElement)
    {
    }
    ElementCollection(const ElementInterface& coordinateElement, DofContainer<const ElementInterface*> dofElement)
        : mCoordinateElement(coordinateElement)
        , mDofElements(dofElement)
    {
    }

    void AddDofElement(const DofType& dofType, const ElementInterface& dofElement)
    {
        mDofElements[dofType] = &dofElement;
    }

    const ElementInterface& CoordinateElement() const
    {
        return mCoordinateElement;
    }

    const ElementInterface& DofElement(const DofType& dofType) const
    {
        return *mDofElements[dofType];
    }

private:
    std::reference_wrapper<const ElementInterface> mCoordinateElement;
    DofContainer<const ElementInterface*> mDofElements;
};
} /* NuTo */
