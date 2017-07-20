#include "mechanics/loads/LoadNodeGroup.h"

NuTo::LoadNodeGroup::LoadNodeGroup(const Group<NodeBase>* rGroup)
    : LoadBase()
    , mGroup(rGroup)
{
}
