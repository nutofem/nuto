#include "mechanics/groups/GroupBase.h"
#include "base/Exception.h"

//! @brief Constructor
NuTo::GroupBase::GroupBase()
{
}

//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
//! @param rNodePtr Node pointer
void NuTo::GroupBase::AddMember(int, NodeBase*)
{
    throw Exception(__PRETTY_FUNCTION__, "Adding a node is only allowed for node groups.");
}

//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
//! @param rNodePtr Node pointer
void NuTo::GroupBase::AddMember(int, ElementBase*)
{
    throw Exception(__PRETTY_FUNCTION__, "Adding an element is only allowed for element groups.");
}

//! @brief replaces a ptr by another one
//! @param rOldPtr
//! @param rNewPtr
void NuTo::GroupBase::ExchangePtr(int, ElementBase*, ElementBase*)
{
    throw Exception(__PRETTY_FUNCTION__, "Exchanging element ptr is only allowed for element groups.");
}

//! @brief replaces a ptr by another one
//! @param rOldPtr
//! @param rNewPtr
void NuTo::GroupBase::ExchangePtr(int, NodeBase*, NodeBase*)
{
    throw Exception(__PRETTY_FUNCTION__, "Exchanging node ptr is only allowed for node groups.");
}
