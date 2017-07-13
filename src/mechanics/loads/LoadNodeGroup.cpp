#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/loads/LoadNodeGroup.h"

NuTo::LoadNodeGroup::LoadNodeGroup(const Group<NodeBase>* rGroup)
    : LoadBase()
    , mGroup(rGroup)
{
}
