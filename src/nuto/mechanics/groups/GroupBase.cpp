#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"

//! @brief Constructor
NuTo::GroupBase::GroupBase()
{
}

//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
//! @param rNodePtr Node pointer
void NuTo::GroupBase::AddMember(NodeBase* rNodePtr)
{
	throw MechanicsException("[NuTo::GroupBase::AddMember] adding a node is only allowed for node groups.");
}

//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
//! @param rNodePtr Node pointer
void NuTo::GroupBase::AddMember(ElementBase* rElementPtr)
{
	throw MechanicsException("[NuTo::GroupBase::AddMember] adding an element is only allowed for element groups.");
}

//! @brief Removes a node from the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
//! @param rNodePtr Node pointer
void NuTo::GroupBase::RemoveMember(NodeBase* rNodePtr)
{
	throw MechanicsException("[NuTo::GroupBase::RemoveMember] removing a node is only allowed for node groups.");
}

//! @brief Removes an element from the group, is only implemented for Element groups in group.cpp, otherwise throws an exception
//! @param rElementPtr Element pointer
void NuTo::GroupBase::RemoveMember(ElementBase* rElementPtr)
{
	throw MechanicsException("[NuTo::GroupBase::RemoveMember] removing an element is only allowed for element groups.");
}

//! @brief check if a group contains the entry
//! @param rElementPtr Element pointer
//! @return TRUE if rMember is in the group, FALSE otherwise
bool NuTo::GroupBase::Contain(NodeBase* rNodePtr)
{
	throw MechanicsException("[NuTo::GroupBase::Contain] looking for a node is only allowed for node groups.");
}

//! @brief check if a group contains the entry
//! @param rElementPtr Element pointer
//! @return TRUE if rMember is in the group, FALSE otherwise
bool NuTo::GroupBase::Contain(ElementBase* rElementPtr)
{
	throw MechanicsException("[NuTo::GroupBase::Contain] looking for an element is only allowed for element groups.");
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
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize GroupBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::GroupBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::GroupBase)
#endif // ENABLE_SERIALIZATION
