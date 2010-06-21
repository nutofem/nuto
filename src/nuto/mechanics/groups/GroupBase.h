// $Id: $
#ifndef GroupBase_H
#define GroupBase_H

#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/groups/GroupEnum.h"

namespace NuTo
{
class StructureBase;
class ElementBase;
class NodeBase;
template <class T>
class Group;

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

#ifdef ENABLE_SERIALIZATION
	//! @brief serializes the class
	//! @param ar         archive
	//! @param version    version
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
	}
#endif // ENABLE_SERIALIZATION

	//! @brief gives the number of group members
	//! @return number of group members
	virtual int GetNumMembers()const=0;

	//! @brief gives the group type
	//! @return group type
	virtual Groups::eGroupId GetType()const=0;

	//! @brief gives the group type as a string
	//! @return group type
	virtual std::string GetTypeString()const=0;

	//! @brief gives the group type
	//! @return group type
	virtual void Info(int rVerboseLevel, const StructureBase* rStructure)const=0;

// *********************************************************************
// * AddMember routine is to be implemented for all new group entities *
// *********************************************************************
	//! @brief Adds a node to the group, is only implemented for Node groups in group.cpp, otherwise throws an exception
	//! @param rNodePtr Node pointer
	virtual void AddMember(NodeBase* rNodePtr);

	//! @brief Adds an element to the group, is only implemented for Element groups, otherwise throws an exception
	//! @param rNodePtr Node pointer
	virtual void AddMember(ElementBase* rNodePtr);

	//! @brief Removes a node from the group
	//! @param rNodePtr Node pointer
	virtual void RemoveMember(NodeBase* rNodePtr);

	//! @brief Removes an element from the group
	//! @param rElementPtr Element pointer
	virtual void RemoveMember(ElementBase* rElementPtr);

	//! @brief check if a group contains the entry
	//! @param rNodePtr Node pointer
    //! @return TRUE if rMember is in the group, FALSE otherwise
	virtual bool Contain(NodeBase* rNodePtr);

	//! @brief check if a group contains the entry
	//! @param rElementPtr Element pointer
    //! @return TRUE if rMember is in the group, FALSE otherwise
	virtual bool Contain(ElementBase* rElementPtr);

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

}; //class definition
}
#endif //GroupBase_H
