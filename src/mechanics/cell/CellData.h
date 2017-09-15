#pragma once

#include "mechanics/cell/PDE_Element.h"


namespace NuTo
{
//! @brief Extracts 'cell data' like nodal values from the cell
//! @remark The life time of objects of this class is limited to the evaluation of one NuTo::Cell. This class will
//! calculate (and cache) all the information that are required once per cell. Classic example: NodeValues. Even if
//! `GetNodeValues()` is called 42 times, they should only be extracted once. Ehrm [TODO]!
template <int TDim>
class CellData
{
public:
    CellData(const PDE_Element<TDim>& element)
        : mElement(element)
    {
    }

    NodeValues ExtractNodeValues(const DofType& dofType) const
    {
        return mElement.ExtractNodeValues(dofType);
    }

private:
    const PDE_Element<TDim>& mElement;
};
} /* NuTo */
