// $Id$
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/loads/LoadNodeGroup.h"


//! @brief constructor
NuTo::LoadNodeGroup::LoadNodeGroup(int rLoadCase, const Group<NodeBase>* rGroup) :
        LoadBase(rLoadCase), mGroup(rGroup)
{

}
