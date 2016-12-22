// $Id$

#include "mechanics/groups/GroupBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/nodes/NodeBase.h"

//! @brief Constructor
NuTo::GroupBase::GroupBase()
{
}

//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
//! @param rNodePtr Node pointer
void NuTo::GroupBase::AddMember(int rId, NodeBase* rNodePtr)
{
	throw MechanicsException("[NuTo::GroupBase::AddMember] adding a node is only allowed for node groups.");
}

//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
//! @param rNodePtr Node pointer
void NuTo::GroupBase::AddMember(int rId, ElementBase* rElementPtr)
{
	throw MechanicsException("[NuTo::GroupBase::AddMember] adding an element is only allowed for element groups.");
}

//! @brief replaces a ptr by another one
//! @param rOldPtr
//! @param rNewPtr
void NuTo::GroupBase::ExchangePtr(int rId, ElementBase* rOldElementPtr, ElementBase* rNewElementPtr)
{
    throw MechanicsException("[NuTo::GroupBase::ExchangePtr] Exchanging element ptr is only allowed for element groups.");
}

//! @brief replaces a ptr by another one
//! @param rOldPtr
//! @param rNewPtr
void NuTo::GroupBase::ExchangePtr(int rId, NodeBase* rOldNodePtr, NodeBase* rNewNodePtr)
{
    throw MechanicsException("[NuTo::GroupBase::ExchangePtr] Exchanging node ptr is only allowed for node groups.");
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::GroupBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::GroupBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::GroupBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::GroupBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::GroupBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::GroupBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::GroupBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize GroupBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(mName);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GroupBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::GroupBase)
#endif // ENABLE_SERIALIZATION
