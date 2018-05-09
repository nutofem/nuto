#pragma once

#include "nuto/mechanics/dofs/DofContainer.h"

namespace NuTo
{
struct DofInfo
{
    DofContainer<int> numIndependentDofs;
    DofContainer<int> numDependentDofs;

    DofContainer<int> numAllDofs()
    {
        DofContainer<int> allDofCount;
        for (auto dof : numIndependentDofs.DofTypes())
            allDofCount[dof] = numIndependentDofs[dof] + numDependentDofs[dof];

        return allDofCount;
    }

    void Merge(DofType dof, DofInfo other)
    {
        numIndependentDofs[dof] = other.numIndependentDofs[dof];
        numDependentDofs[dof] = other.numDependentDofs[dof];
    }
};
} /* NuTo */
