#pragma once

#include "mechanics/elements/ElementCollection.h"
#include "mechanics/nodes/DofContainer.h"

namespace NuTo
{


//! @brief Extracts 'cell data' like nodal values from the cell
//! @remark The life time of objects of this class is limited to the evaluation of one NuTo::Cell. This class will
//! calculate (and cache) all the information that are required once per cell. Classic example: NodeValues. Even if
//! `GetNodeValues()` is called 42 times, they should only be extracted once.
class CellData
{
public:
    CellData(const ElementCollection& elements)
        : mElements(elements)
    {
    }

    NodeValues GetCoordinates() const
    {
        return mElements.CoordinateElement().ExtractNodeValues();
    }

    NodeValues GetNodeValues(const DofType& dofType) const
    {
        NodeValues& nodeValues = mNodeValues[dofType];
        if (nodeValues.size() == 0)
        {
            /* The mutable member mNodeValues may be accessed by multiple threads
             * at once - in an OpenMP context.
             * Concurrent read access is no problem. This if branch, however,
             * performs write access. No concurrency here.
             */
#ifdef _OPENMP
#pragma omp critical
#endif
            nodeValues = mElements.DofElement(dofType).ExtractNodeValues();
        }
        return nodeValues;
    }

private:
    mutable DofContainer<NodeValues> mNodeValues;
    const ElementCollection& mElements;
};
} /* NuTo */
