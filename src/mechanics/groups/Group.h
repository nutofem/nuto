#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/map.hpp>
#include <boost/utility/identity_type.hpp>
#else
#include <map>
#endif // ENABLE_SERIALIZATION

#include "mechanics/groups/GroupBase.h"
#include "mechanics/MechanicsException.h"

namespace NuTo
{
class ElementBase;
class NodeBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all groups
template <class T>
class Group : public GroupBase, public std::map<int,T*>
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    Group() : GroupBase(), std::map<int, T*>(){}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes (saves) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start saving Group<T>" << std::endl;
#endif
        // save the Addresses in another map
        std::map<int, std::uintptr_t> mapCast;
        for (typename std::map<int,T*>::const_iterator it = this->begin(); it!= this->end(); it++)
        {
            mapCast.insert(std::pair<int, std::uintptr_t>(it->first, reinterpret_cast<std::uintptr_t>(it->second)));
        }

        // serialize this map containing Addresses
        ar & boost::serialization::make_nvp ("map", mapCast );
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GroupBase);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish saving Group<T>" << std::endl;
#endif
    }

    //! @brief deserializes (loads) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start loading Group<T>" << std::endl;
#endif
        std::map<int, std::uintptr_t> mapCast;
        ar & boost::serialization::make_nvp ("map", mapCast );
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GroupBase);

        for (std::map<int,std::uintptr_t>::iterator it = mapCast.begin(); it!= mapCast.end(); it++)
        {
            this->insert(std::pair<int, T*>(it->first, reinterpret_cast<T*>(it->second)));
        }
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish loading Group<T>" << std::endl;
#endif
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

//    //! @brief deserializes (loads) the class
//    //! @param ar         archive
//    //! @param version    version
//    template<class Archive>
//    void serialize(Archive & ar, const unsigned int version)
//    {
//#ifdef DEBUG_SERIALIZATION
//        std::cout << "start loading Group<T>" << std::endl;
//#endif
//        std::map<int,std::uintptr_t>* mapCast = dynamic_cast<std::map<int,std::uintptr_t>* >(this);
//        ar & boost::serialization::make_nvp ("map", *mapCast );
//        dynamic_cast<std::map<int,T*>* >(this) = mapCast;

//        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GroupBase);


//#ifdef DEBUG_SERIALIZATION
//        std::cout << "finish loading Group<T>" << std::endl;
//#endif
//    }

    void SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast) override
    {
        for (typename std::map<int,T*>::iterator it = this->begin(); it!= this->end(); it++)
        {
            std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast =
                    mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(it->second));
            if(itCast != mNodeMapCast.end())
            {
                it->second = reinterpret_cast<T*>(itCast->second);
            }
            else
                throw MechanicsException("[NuTo::Group] The NodeBase/ElementBase-Pointer could not be updated.");
        }
    }
#endif // ENABLE_SERIALIZATION
    
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
    	if ((this->insert(std::pair<int,T*>(rId,rMember))).second==false)
            throw MechanicsException("[Group::AddMember] Group member already exists in the group.");
    }

    //! @brief removes a group member
    //! @param rMember member to be removed
    void RemoveMember(int rId) override
    {
        if (!this->erase(rId))
            throw MechanicsException("[Group::AddMember] Group member to be deleted is not within the group.");
    }

    //! @brief check if a group contains the entry
    //! @param rElementPtr Element pointer
    //! @return TRUE if rMember is in the group, FALSE otherwise
    bool Contain(int rId) const override
    {
        if (this->find(rId)==this->end())
        	return false;
        else
        	return true;
    }


    //! @brief replaces a ptr by another one
    //! @param rOldPtr
    //! @param rNewPtr
    void ExchangePtr(int rId, T* rOldMember, T* rNewMember) override
    {
        typename std::map<int,T*>::iterator it(this->find(rId));
        if (it==this->end())
        	throw MechanicsException("[Group::ExchangePtr] New group member can not be inserted, since the id does not exist in group.");
        else
        {
        	assert(it->second==rOldMember);
        	it->second=rNewMember;
        }
    }

    //! @brief joins two groups
    //! @return group
    GroupBase* Unite (const NuTo::GroupBase* rOther) const override
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==nullptr)
			throw MechanicsException("[NuTo::Group::Unite] Groups do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_union (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

    //! @brief returns a group with all members of current group, which are not presented in the second group
    //! @return group
    GroupBase* Difference (const NuTo::GroupBase* rOther) const override
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==nullptr)
			throw MechanicsException("[NuTo::Group::Difference] Groups do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_difference (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

    //! @brief returns a group with all members which are elements of both groups
    //! @return group
    GroupBase* Intersection (const NuTo::GroupBase* rOther) const override
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==nullptr)
			throw MechanicsException("[NuTo::Group::Intersection] Groups do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_intersection (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

    //! @brief returns a group with the symmetric difference, i.e. all elements present in group 1 or in group 2 (but not in both)
    //! @return group
    GroupBase* SymmetricDifference (const NuTo::GroupBase* rOther) const override
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==nullptr)
			throw MechanicsException("[NuTo::Group::SymmetricDifference] Groups to be united do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_symmetric_difference (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

    //! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
    Group<ElementBase>* AsGroupElement() override;

    //! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
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
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::Group<NuTo::NodeBase>)
BOOST_CLASS_EXPORT_KEY(NuTo::Group<NuTo::ElementBase>)
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((std::map<int,NuTo::NodeBase*>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((std::map<int,NuTo::ElementBase*>)))
#endif // SWIG
#endif // ENABLE_SERIALIZATION

