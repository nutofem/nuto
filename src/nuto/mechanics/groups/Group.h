#ifndef GROUP_H
#define GROUP_H
#include <algorithm>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/set.hpp>
#else
#include <set>
#endif // ENABLE_SERIALIZATION


#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class ElementBase;
class NodeBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all groups
template <class T>
class Group : public GroupBase, public std::set<T*>
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    Group() : GroupBase(), std::set<T*>(){};

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize Group<T>" << std::endl;
#endif
        ar & boost::serialization::make_nvp ("set",boost::serialization::base_object< std::set<T*> > ( *this ) )
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GroupBase);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize Group<T>" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION
    
    //! @brief gives the number of group members
    //! @return number of group members
    int GetNumMembers()const
    {
        return (int)this->size();
    }

    //! @brief adds a group member
    //! @param rMember new member
    void AddMember(T* rMember)
    {
        if (!this->insert(rMember).second)
            throw MechanicsException("[Group::AddMember] Group member already exists in the group.");
    }

    //! @brief removes a group member
    //! @param rMember member to be removed
    void RemoveMember(T* rMember)
    {
        if (!this->erase(rMember))
            throw MechanicsException("[Group::AddMember] Group member to be deleted is not within the group.");
    }
    
	//! @brief check if a group contains the entry
	//! @param rElementPtr Element pointer
    //! @return TRUE if rMember is in the group, FALSE otherwise
    bool Contain(T* rMember)
    {
        if (this->find(rMember)==this->end())	return false;
        else							       	return true;
    }

    //! @brief replaces a ptr by another one
    //! @param rOldPtr
    //! @param rNewPtr
    void ExchangePtr(T* rOldMember, T* rNewMember)
    {
        if (this->erase(rOldMember))
        {
            if (this->insert(rNewMember).second==false)
            {
                throw MechanicsException("[Group::ExchangePtr] New group member can not be inserted.");
            }
        }
    }

    //! @brief joins two groups
    //! @return group
    GroupBase* Unite (const NuTo::GroupBase* rOther)const
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==0)
			throw MechanicsException("[NuTo::Group::Unite] Groups do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_union (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

    //! @brief returns a group with all members of current group, which are not presented in the second group
    //! @return group
    GroupBase* Difference (const NuTo::GroupBase* rOther)const
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==0)
			throw MechanicsException("[NuTo::Group::Difference] Groups do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_difference (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

    //! @brief returns a group with all members which are elements of both groups
    //! @return group
    GroupBase* Intersection (const NuTo::GroupBase* rOther)const
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==0)
			throw MechanicsException("[NuTo::Group::Intersection] Groups do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_intersection (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

    //! @brief returns a group with the symmetric difference, i.e. all elements present in group 1 or in group 2 (but not in both)
    //! @return group
    GroupBase* SymmetricDifference (const NuTo::GroupBase* rOther)const
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==0)
			throw MechanicsException("[NuTo::Group::SymmetricDifference] Groups to be united do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_symmetric_difference (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

    //! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
    Group<ElementBase>* AsGroupElement();

    //! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
    const Group<ElementBase>* AsGroupElement()const;

    //! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
    Group<NodeBase>* AsGroupNode();

    //! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
    const Group<NodeBase>* AsGroupNode()const;

    //! @brief gives the group type
    //! @return group type
    Groups::eGroupId GetType()const;

    //! @brief gives the group type
    //! @return group type as string
    std::string GetTypeString()const;

    //! @brief gives the group type
    //! @return group type
    void Info(int rVerboseLevel, const NuTo::StructureBase* rStructure)const;

};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::Group<NuTo::NodeBase>)
BOOST_CLASS_EXPORT_KEY(NuTo::Group<NuTo::ElementBase>)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
#endif //GROUP_H

