// $Id$
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/loads/LoadNode.h"


//! @brief constructor
NuTo::LoadNode::LoadNode(int rLoadCase, const NodeBase* rNode) : LoadBase(rLoadCase), mNode(rNode)
{}
