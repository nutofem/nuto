#pragma once

#include "nuto/base/Group.h"

#include "nuto/mechanics/dofs/DofInfo.h"
#include "nuto/mechanics/constraints/Constraints.h"
#include "nuto/mechanics/nodes/DofNode.h"

namespace NuTo
{

namespace DofNumbering
{

//! Builds the dof numbering for type `dof` for all the `dofNodes` and stores it at the `dofNodes`.
//! @remark Numbering starts at 0 with all the dofs that are not constrained. Constrained dofs have the highest numbers
DofInfo Build(const Group<DofNode>& dofNodes, DofType dof, const Constraint::Constraints& constraints);

//! Extracts the dof numbers from `dofNodes` for a given component (0 for scalar, x=0, y=1, z=2 for vector dofs)
std::vector<int> Get(const Group<DofNode>& dofNodes, int component = 0);

} /* DofNumbering */
} /* NuTo */
