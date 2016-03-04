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
#include "nuto/mechanics/nodes/NodeDof.h"
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
// serializes the class
template void NuTo::ConstraintNode::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNode::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNode::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNode::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNode::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNode::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintNode::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintNode" << std::endl;
#endif
    std::uintptr_t& mNodeAdress = reinterpret_cast<std::uintptr_t&>(mNode);
    ar & boost::serialization::make_nvp("mNode", mNodeAdress);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintNode" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintNode)

void NuTo::ConstraintNode::SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
{
    int t = reinterpret_cast<std::uintptr_t>(mNode);
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
