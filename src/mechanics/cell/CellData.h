#pragma once

#include "mechanics/elements/ElementSimple.h"
#include "mechanics/nodes/DofContainer.h"

namespace NuTo
{


//! @brief Extracts 'cell data' like nodal values from the cell
//! @remark The life time of objects of this class is limited to the evaluation of one NuTo::Cell. This class will
//! calculate (and cache) all the information that are required once per cell. Classic example: NodeValues. Even if
//! `GetNodeValues()` is called 42 times, they should only be extracted once. Ehrm [TODO]!
class CellData
{
public:
    CellData(const DofContainer<ElementSimple*>& rElements)
        : mElements(rElements)
    {
    }

    NodeValues GetNodeValues(const DofType& rDofType) const
    {
        return mElements[rDofType]->ExtractNodeValues();
    }

private:
    const DofContainer<ElementSimple*>& mElements;
};
} /* NuTo */
