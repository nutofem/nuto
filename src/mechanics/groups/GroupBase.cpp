
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
    throw Exception("[NuTo::GroupBase::AddMember] adding a node is only allowed for node groups.");
}

//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
//! @param rNodePtr Node pointer
void NuTo::GroupBase::AddMember(int rId, ElementBase* rElementPtr)
{
    throw Exception("[NuTo::GroupBase::AddMember] adding an element is only allowed for element groups.");
}

//! @brief replaces a ptr by another one
//! @param rOldPtr
//! @param rNewPtr
void NuTo::GroupBase::ExchangePtr(int rId, ElementBase* rOldElementPtr, ElementBase* rNewElementPtr)
{
    throw Exception(
            "[NuTo::GroupBase::ExchangePtr] Exchanging element ptr is only allowed for element groups.");
}

//! @brief replaces a ptr by another one
//! @param rOldPtr
//! @param rNewPtr
void NuTo::GroupBase::ExchangePtr(int rId, NodeBase* rOldNodePtr, NodeBase* rNewNodePtr)
{
    throw Exception("[NuTo::GroupBase::ExchangePtr] Exchanging node ptr is only allowed for node groups.");
}
