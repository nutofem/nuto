#ifndef GROUP_H
#define GROUP_H
#include <set>
#include <algorithm>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
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
//        ar & BOOST_SERIALIZATION_NVP(mMembers);
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
    
    //! @brief joins two groups
    //! @return group
    GroupBase* Unite (const NuTo::GroupBase* rOther)const
    {
        NuTo::Group<T>* returnGroup = new NuTo::Group<T>();
		const Group<T>* rOtherT = dynamic_cast<const Group<T>* >(rOther);
		if (rOtherT==0)
			throw MechanicsException("[NuTo::Group::Unite] Groups to be united do not have the same type.");
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
			throw MechanicsException("[NuTo::Group::Unite] Groups to be united do not have the same type.");
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
			throw MechanicsException("[NuTo::Group::Unite] Groups to be united do not have the same type.");
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
			throw MechanicsException("[NuTo::Group::Unite] Groups to be united do not have the same type.");
		std::insert_iterator<NuTo::Group<T> > returnGroupInsertIterator (*returnGroup, returnGroup->begin ());
		std::set_symmetric_difference (this->begin(), this->end(), rOtherT->begin(), rOtherT->end(), returnGroupInsertIterator);
		return returnGroup;
    }

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
#endif //GROUP_H

