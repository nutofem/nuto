#pragma once

#include "base/UniqueId.h"
#include "mechanics/dofs/DofContainer.h"
#include "mechanics/elements/ElementCollection.h"

namespace NuTo
{


//! @brief Extracts 'cell data' like nodal values from the cell
//! @remark The life time of objects of this class is limited to the evaluation of one NuTo::Cell. This class will
//! calculate (and cache) all the information that are required once per cell. Classic example: NodeValues. Even if
//! `GetNodeValues()` is called 42 times, they should only be extracted once.
class CellData : public UniqueId<CellData>
{
public:
    CellData(const ElementCollection& elements, int numIntegrationPoints)
        : mElements(elements)
        , mNumIntegrationPoints(numIntegrationPoints)
    {
    }

    int GetNumIntegrationPoints() const
    {
        return mNumIntegrationPoints;
    }

    NodeValues GetCoordinates() const
    {
        return mElements.CoordinateElement().ExtractNodeValues();
    }

    NodeValues GetNodeValues(const DofType& dofType) const
    {
        NodeValues& nodeValues = mNodeValues[dofType];
        if (nodeValues.size() == 0)
            nodeValues = mElements.DofElement(dofType).ExtractNodeValues();

        return nodeValues;
    }

private:
    mutable DofContainer<NodeValues> mNodeValues;
    const ElementCollection& mElements;
    int mNumIntegrationPoints;
};
} /* NuTo */
