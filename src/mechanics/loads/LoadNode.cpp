#include "mechanics/nodes/NodeBase.h"
#include "mechanics/loads/LoadNode.h"

NuTo::LoadNode::LoadNode(const NodeBase* rNode)
    : LoadBase()
    , mNode(rNode)
{
}
