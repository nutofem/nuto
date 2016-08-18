// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constraints/ConstraintNode.h"


//! @brief constructor
NuTo::ConstraintNode::ConstraintNode(const NodeBase* rNode) : mNode(rNode)
{
}

//! @brief destructor
NuTo::ConstraintNode::~ConstraintNode()
{
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintNode)
void NuTo::ConstraintNode::SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
{
    std::map<std::uintptr_t, std::uintptr_t>::const_iterator it = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(mNode));
    if (it!=mNodeMapCast.end())
    {
        NodeBase** temp = const_cast<NodeBase**>(&mNode);
        *temp = reinterpret_cast<NodeBase*>(it->second);
    }
    else
        throw MechanicsException("[NuTo::ConstraintNode] The NodeBase-Pointer could not be updated.");
}

#endif // ENABLE_SERIALIZATION
