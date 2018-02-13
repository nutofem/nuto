#pragma once

#include "base/Group.h"

#include "mechanics/dofs/DofInfo.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/nodes/NodeSimple.h"

namespace NuTo
{

namespace DofNumbering
{

DofInfo Build(const Group<NodeSimple>& dofNodes, DofType dof, const Constraint::Constraints& constraints);

} /* DofNumbering */
} /* NuTo */
