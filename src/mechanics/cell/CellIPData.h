#pragma once

#include "mechanics/elements/Element.h"
#include "mechanics/nodes/DofContainer.h"
#include "mechanics/cell/Jacobian.h"

namespace NuTo
{


//! @brief Similar to NuTo::CellData, but for N and B
class CellIPData
{
public:
    CellIPData(const DofContainer<const Element*> rElements, const NuTo::Jacobian& rJacobian)
        : mElements(rElements)
        , mJacobian(rJacobian)
    {
    }


private:
    const DofContainer<const Element*> mElements;
    const NuTo::Jacobian& mJacobian;
};
};
} /* NuTo */
