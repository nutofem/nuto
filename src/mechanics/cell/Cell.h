#pragma once

#include <vector>
#include "mechanics/elements/Element.h"
#include "mechanics/nodes/DofContainer.h"

namespace NuTo
{
class Cell
{
public:
    Cell(const Element& rCoordinateElement)
        : mCoordinateElement(rCoordinateElement)
    {
    }

private:
    const Element& mCoordinateElement;
    DofContainer<Element*> mElements;
}; 
} /* NuTo */ 
