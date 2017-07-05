#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include <cassert>

using namespace NuTo;

int NodeBase::GetDof(Node::eDof rDof) const
{
    assert(GetNum(rDof) == 1);
    return GetDof(rDof, 0);
}

namespace NuTo
{
std::ostream& operator<<(std::ostream& out, const NodeBase& node)
{
    node.Info(out);
    return out;
}
} /* NuTo */

