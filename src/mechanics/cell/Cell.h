#pragma once

#include <vector>
#include "mechanics/elements/Element.h"

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
}; 
} /* NuTo */ 
