#pragma once

#include "base/Group.h"

#include "mechanics/constraints/Constraints.h"
#include "mechanics/dofs/DofContainer.h"
#include "mechanics/nodes/NodeSimple.h"

namespace NuTo
{

namespace DofNumbering
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

DofInfo Build(const Groups::Group<NodeSimple>& dofNodes, DofType dof, const Constraint::Constraints& constraints);

} /* DofNumbering */
} /* NuTo */
