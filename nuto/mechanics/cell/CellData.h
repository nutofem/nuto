#pragma once

#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/elements/ElementCollection.h"

namespace NuTo
{


//! @brief Extracts 'cell data' like nodal values from the cell
//! @remark The life time of objects of this class is limited to the evaluation of one NuTo::Cell. This class will
//! calculate (and cache) all the information that are required once per cell. Classic example: NodeValues. Even if
//! `GetNodeValues()` is called 42 times, they should only be extracted once.
class CellData
{
public:
    CellData(const ElementCollection& elements, int cellId)
        : mElements(elements)
        , mCellId(cellId)
    {
    }

    int GetCellId() const
    {
        return mCellId;
    }

    NodeValues GetCoordinates() const
    {
        return mElements.CoordinateElement().ExtractNodeValues();
    }

    const NodeValues& GetNodeValues(DofType dofType) const
    {
        NodeValues& nodeValues = mNodeValues[dofType];
        if (nodeValues.size() == 0)
            nodeValues = mElements.DofElement(dofType).ExtractNodeValues();

        return nodeValues;
    }

    const ElementCollection& Elements() const
    {
        return mElements;
    }

private:
    mutable DofContainer<NodeValues> mNodeValues;
    const ElementCollection& mElements;
    int mCellId;
};
} /* NuTo */
