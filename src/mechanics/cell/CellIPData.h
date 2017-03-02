#pragma once

#include "mechanics/elements/Element.h"
#include "mechanics/nodes/DofContainer.h"

namespace NuTo
{


//! @brief Extracts 'cell ip data' like N and B matrix from the cell
class CellIPData
{
    public:
    CellIPData(const DofContainer<const Element*> rElements)
        : mElements(rElements)
    {
    }



private:
    const DofContainer<const Element*> mElements;
};
};
} /* NuTo */
