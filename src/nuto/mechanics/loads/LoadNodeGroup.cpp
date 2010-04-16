// $Id$
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/loads/LoadNodeGroup.h"


//! @brief constructor
NuTo::LoadNodeGroup::LoadNodeGroup(const Group<NodeBase>* rGroup) :
        LoadBase(), mGroup(rGroup)
{

}
