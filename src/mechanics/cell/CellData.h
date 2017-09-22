#pragma once

#include "mechanics/interpolation/CellInterpolationBase.h"
#include "mechanics/dofs/DofContainer.h"

namespace NuTo
{


//! @brief Extracts 'cell data' like nodal values from the cell
//! @remark The life time of objects of this class is limited to the evaluation of one NuTo::Cell. This class will
//! calculate (and cache) all the information that are required once per cell. Classic example: NodeValues. Even if
//! `GetNodeValues()` is called 42 times, they should only be extracted once.
class CellData
{
public:
    CellData(const DofContainer<CellInterpolationBase*>& elements)
        : mElements(elements)
    {
    }

    NodeValues GetNodeValues(const DofType& dofType) const
    {
        NodeValues& nodeValues = mNodeValues[dofType];
        if (nodeValues.size() == 0)
        {
#ifdef _OPENMP
#pragma omp critical
#endif
            nodeValues = mElements[dofType]->ExtractNodeValues();
        }
        return nodeValues;
    }

private:
    mutable DofContainer<NodeValues> mNodeValues;
    const DofContainer<CellInterpolationBase*>& mElements;
};
} /* NuTo */
