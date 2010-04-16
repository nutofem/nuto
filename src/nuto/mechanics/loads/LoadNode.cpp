// $Id$
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/loads/LoadNode.h"


//! @brief constructor
NuTo::LoadNode::LoadNode(const NodeBase* rNode) : LoadBase(), mNode(rNode)
{}
