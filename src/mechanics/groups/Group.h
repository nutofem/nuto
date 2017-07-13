#pragma once

#include <map>

#include <algorithm>
#include "mechanics/groups/GroupBase.h"
#include "base/Exception.h"

namespace NuTo
{
class ElementBase;
class NodeBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all groups
template <class T>
class Group : public GroupBase, public std::map<int, T*>
{

public:
    //! @brief constructor
    Group()
        : GroupBase()
        , std::map<int, T*>()
    {
    }

    //! @brief gives the number of group members
    //! @return number of group members
    int GetNumMembers() const override
    {
        return (int)this->size();
    }

    //! @brief returns the group members
    //! @return group members (id)
    std::vector<int> GetMemberIds() const override
    {
        std::vector<int> members;
        members.reserve(this->size());
        for (auto pair : *this)
            members.push_back(pair.first);
        return members;
    }

    //! @brief adds a group member
    //! @param rMember new member
    void AddMember(int rId, T* rMember) override
    {
        if ((this->insert(std::pair<int, T*>(rId, rMember))).second == false)
            throw Exception("[Group::AddMember] Group member already exists in the group.");
    }

    //! @brief removes a group member
    //! @param rMember member to be removed
    void RemoveMember(int rId) override
    {
        if (!this->erase(rId))
            throw Exception("[Group::AddMember] Group member to be deleted is not within the group.");
    }

    //! @brief check if a group contains the entry
    //! @param rElementPtr Element pointer
    //! @return TRUE if rMember is in the group, FALSE otherwise
    bool Contain(int rId) const override
    {
        return this->find(rId) != this->end();
    }


    //! @brief replaces a ptr by another one
    //! @param rOldPtr
    //! @param rNewPtr
    void ExchangePtr(int rId, T* rOldMember, T* rNewMember) override
    {
        typename std::map<int, T*>::iterator it(this->find(rId));
        if (it == this->end())
            throw Exception(
                    "[Group::ExchangePtr] New group member can not be inserted, since the id does not exist in group.");
        else
        {
            assert(it->second == rOldMember);
            it->second = rNewMember;
        }
    }

    //! @brief joins two groups
    //! @return group
    GroupBase* Unite(const NuTo::GroupBase* rOther) const override
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
        const Group<T>* rOtherT = dynamic_cast<const Group<T>*>(rOther);
        if (rOtherT == nullptr)
            throw Exception("[NuTo::Group::Unite] Groups do not have the same type.");
        std::insert_iterator<NuTo::Group<T>> returnGroupInsertIterator(*returnGroup, returnGroup->begin());
        std::set_union(this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
        return returnGroup;
    }

    //! @brief Union of groupOne and groupTwo
    //! @return Union of groupOne and groupTwo
    static Group<T> Unite(const Group<T>& groupOne, const Group<T>& groupTwo)
    {
        NuTo::Group<T> returnGroup;
        std::insert_iterator<NuTo::Group<T>> returnGroupInsertIterator(returnGroup, returnGroup.begin());
        std::set_union(groupOne.begin(), groupOne.end(), groupTwo.begin(), groupTwo.end(), returnGroupInsertIterator);
        return returnGroup;
    }

    //! @brief returns a group with all members of current group, which are not presented in the second group
    //! @return group
    GroupBase* Difference(const NuTo::GroupBase* rOther) const override
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
        const Group<T>* rOtherT = dynamic_cast<const Group<T>*>(rOther);
        if (rOtherT == nullptr)
            throw Exception("[NuTo::Group::Difference] Groups do not have the same type.");
        std::insert_iterator<NuTo::Group<T>> returnGroupInsertIterator(*returnGroup, returnGroup->begin());
        std::set_difference(this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
        return returnGroup;
    }

    //! @brief returns a group with all members which are elements of both groups
    //! @return group
    GroupBase* Intersection(const NuTo::GroupBase* rOther) const override
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
        const Group<T>* rOtherT = dynamic_cast<const Group<T>*>(rOther);
        if (rOtherT == nullptr)
            throw Exception("[NuTo::Group::Intersection] Groups do not have the same type.");
        std::insert_iterator<NuTo::Group<T>> returnGroupInsertIterator(*returnGroup, returnGroup->begin());
        std::set_intersection(this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
        return returnGroup;
    }

    //! @brief returns a group with the symmetric difference, i.e. all elements present in group 1 or in group 2 (but
    //! not in both)
    //! @return group
    GroupBase* SymmetricDifference(const NuTo::GroupBase* rOther) const override
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
        const Group<T>* rOtherT = dynamic_cast<const Group<T>*>(rOther);
        if (rOtherT == nullptr)
            throw Exception(
                    "[NuTo::Group::SymmetricDifference] Groups to be united do not have the same type.");
        std::insert_iterator<NuTo::Group<T>> returnGroupInsertIterator(*returnGroup, returnGroup->begin());
        std::set_symmetric_difference(this->begin(), this->end(), rOtherT->begin(), rOtherT->end(),
                                      returnGroupInsertIterator);
        return returnGroup;
    }

    //! @brief either casts the pointer to an element group or throws an exception for groups which are not element
    //! groups
    Group<ElementBase>* AsGroupElement() override;

    //! @brief either casts the pointer to an element group or throws an exception for groups which are not element
    //! groups
    const Group<ElementBase>* AsGroupElement() const override;

    //! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
    Group<NodeBase>* AsGroupNode() override;

    //! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
    const Group<NodeBase>* AsGroupNode() const override;

    //! @brief gives the group type
    //! @return group type
    eGroupId GetType() const override;

    //! @brief gives the group type
    //! @return group type as string
    std::string GetTypeString() const override;

    //! @brief gives the group type
    //! @return group type
    void Info(int rVerboseLevel, const NuTo::StructureBase* rStructure) const override;
};
} // namespace NuTo
