#include "mechanics/nodes/NodeBase.h"
#include "mechanics/loads/LoadNode.h"

NuTo::LoadNode::LoadNode(const NodeBase* rNode)
    : LoadBase()
    , mNode(rNode)
{
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadNode)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::LoadNode)
BOOST_CLASS_TRACKING(NuTo::LoadNode, track_always)
#endif
