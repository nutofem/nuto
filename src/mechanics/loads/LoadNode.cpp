// $Id$
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/loads/LoadNode.h"


//! @brief constructor
NuTo::LoadNode::LoadNode(int rLoadCase, const NodeBase* rNode) : LoadBase(rLoadCase), mNode(rNode)
{}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadNode)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::LoadNode)
BOOST_CLASS_TRACKING(NuTo::LoadNode, track_always)
#endif
