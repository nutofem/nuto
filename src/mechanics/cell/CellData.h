#pragma once

#include "mechanics/elements/Element.h"
#include "mechanics/nodes/DofContainer.h"

namespace NuTo
{


//! @brief Extracts 'cell data' like nodal values from the cell
class CellData
{
public:
    CellData(const DofContainer<const Element*> rElements)
        : mElements(rElements)
    {
    }

    Eigen::VectorXd GetNodeValues(const DofType& rDofType)
    {
        return mElements[rDofType]->ExtractNodeValues();
    }

private:
    const DofContainer<const Element*> mElements;
};
} /* NuTo */
