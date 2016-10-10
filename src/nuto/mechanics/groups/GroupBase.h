// $Id$
#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <map>
#endif // ENABLE_SERIALIZATION

#include <eigen3/Eigen/Core>
#include <string>

namespace NuTo
{
class StructureBase;
class ElementBase;
class NodeBase;
template <class T> class Group;
template <class T, int rows> class FullVector;
enum class eGroupId;

//! @author Jörg F. Unger, ISM
//! @date November 2009
//! @brief ... standard abstract class for all groups
class GroupBase
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
	//! @brief constructor
	GroupBase();

	//! @brief ... destructor
    virtual ~GroupBase(){}

#ifdef ENABLE_SERIALIZATION
	//! @brief serializes the class
	//! @param ar         archive
	//! @param version    version
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast   std::map containing the old and new Addresses
    virtual void SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)=0;
#endif // ENABLE_SERIALIZATION

	//! @brief gives the number of group members
	//! @return number of group members
	void SetName(const std::string& rName)
	{
		mName = rName;
	}

	//! @brief gives the number of group members
	//! @return number of group members
	std::string GetName()const
	{
		return mName;
	}

	//! @brief gives the number of group members
	//! @return number of group members
	virtual int GetNumMembers()const=0;

	//! @brief gives the group member ids
	//! @return group members
	virtual NuTo::FullVector<int,Eigen::Dynamic> GetMemberIds()const=0;

	//! @brief gives the group type
	//! @return group type
    virtual eGroupId GetType()const=0;

	//! @brief gives the group type as a string
	//! @return group type
	virtual std::string GetTypeString()const=0;

	//! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
	virtual Group<ElementBase>* AsGroupElement()=0;

    //! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
	virtual const Group<ElementBase>* AsGroupElement()const=0;

    //! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
	virtual Group<NodeBase>* AsGroupNode()=0;

    //! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
	virtual const Group<NodeBase>* AsGroupNode()const=0;

	//! @brief gives the group type
	//! @return group type
	virtual void Info(int rVerboseLevel, const StructureBase* rStructure)const=0;

// *********************************************************************
// * AddMember routine is to be implemented for all new group entities *
// *********************************************************************
	//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
	//! @param rNodePtr Node pointer
	virtual void AddMember(int rId, NodeBase* rNodePtr);

	//! @brief Adds an element to the group, is only implemented for Element groups, otherwise throws an exception
	//! @param rNodePtr Node pointer
	virtual void AddMember(int rId, ElementBase* rNodePtr);

	//! @brief Removes a node from the group
	//! @param rNodePtr Node pointer
	virtual void RemoveMember(int rId)=0;

	//! @brief check if a group contains the entry
	//! @param rId id of the entry
    //! @return TRUE if rMember is in the group, FALSE otherwise
	virtual bool Contain(int rId)const=0;

    //! @brief replaces a ptr by another one
    //! @param rOldPtr
    //! @param rNewPtr
    virtual void ExchangePtr(int rId, ElementBase* rOldElementPtr, ElementBase* rNewElementPtr);

    //! @brief replaces a ptr by another one
    //! @param rOldPtr
    //! @param rNewPtr
    virtual void ExchangePtr(int rId, NodeBase* rOldNodePtr, NodeBase* rNewNodePtr);

	//! @brief Unites two groups
	//! @param rOther other group
	//! @return newly created united group
	virtual GroupBase* Unite (const NuTo::GroupBase* rOther)const=0;

	//! @brief Difference of two groups
	//! @param rOther other group
	//! @return newly created united group
	virtual GroupBase* Difference (const NuTo::GroupBase* rOther)const=0;

	//! @brief Intersection of two groups
	//! @param rOther other group
	//! @return newly created united group
	virtual GroupBase* Intersection (const NuTo::GroupBase* rOther)const=0;

	//! @brief Symmetric difference of two groups
	//! @param rOther other group
	//! @return newly created united group
	virtual GroupBase* SymmetricDifference (const NuTo::GroupBase* rOther)const=0;

protected:
    //! @brief ... name of the group
	std::string mName;

}; //class definition
}
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::GroupBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

