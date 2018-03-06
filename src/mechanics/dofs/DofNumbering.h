#pragma once

#include "base/Group.h"

#include "mechanics/dofs/DofInfo.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/nodes/NodeSimple.h"

namespace NuTo
{

namespace DofNumbering
{

//! Builds the dof numbering for type `dof` for all the `dofNodes` and stores it at the `dofNodes`.
//! @remark Numbering starts at 0 with all the dofs that are not constrained. Constraned dofs have the highest numbers
DofInfo Build(const Group<NodeSimple>& dofNodes, DofType dof, const Constraint::Constraints& constraints);

//! Extracts the dof numbers from `dofNodes` for a given component (0 for scalar, x=0, y=1, z=2 for vector dofs)
std::vector<int> Get(const Group<NodeSimple>& dofNodes, int component = 0);

} /* DofNumbering */
} /* NuTo */
