#pragma once

#include "nuto/mechanics/dofs/DofContainer.h"

namespace NuTo
{
struct DofInfo
{
    DofContainer<int> numIndependentDofs;
    DofContainer<int> numDependentDofs;

    void Merge(DofType dof, DofInfo other)
    {
        numIndependentDofs[dof] = other.numIndependentDofs[dof];
        numDependentDofs[dof] = other.numDependentDofs[dof];
    }
};
} /* NuTo */
