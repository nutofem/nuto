// $Id: $

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"

void NuTo::StructureBase::GroupInfo(int rVerboseLevel)const
{
    std::cout << "number of groups  : " << mGroupMap.size() << std::endl;
    if (rVerboseLevel>2)
    {
        for (boost::ptr_map<std::string,GroupBase>::const_iterator it= mGroupMap.begin(); it!=mGroupMap.end(); it++)
        {
            std::cout << "  Group " << it->first << std::endl;
            it->second->Info(rVerboseLevel, this);
            std::cout <<  std::endl;
        }
    }

}

//! @brief gives the identifier of a group
//! @param pointer to a group
//! @return identifier
std::string NuTo::StructureBase::GroupGetId(GroupBase* rGroup)const
{
    for (boost::ptr_map<std::string,GroupBase>::const_iterator
            it = mGroupMap.begin(); it!= mGroupMap.end(); it++)
    {
        if (it->second==rGroup)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::GroupGetId] Group does not exist.");
}

//! @brief ... Creates a group for the structure
//! @param ... rIdent identifier for the group
//! @param ... rType  type of the group, e.g. "NODES" or "ELEMENTS"
void NuTo::StructureBase::GroupCreate(const std::string& rIdent, const std::string& rType)
{
    // transform string to uppercase
    std::string GroupTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(GroupTypeString), (int(*)(int)) toupper);

    if (GroupTypeString==std::string("ELEMENTS"))
    {
        mGroupMap.insert(const_cast<std::string&>(rIdent), new NuTo::Group<ElementBase>);
    }
    else
    {
		if (GroupTypeString==std::string("NODES"))
		{
			mGroupMap.insert(const_cast<std::string&>(rIdent), new NuTo::Group<NodeBase>);
		}
		else
		{
			throw MechanicsException("[NuTo::StructureBase::GroupCreate] Group type not implemented.");
		}
    }
}

//! @brief ... Adds a node to a node group
//! @param ... rIdentGroup identifier for the group
//! @param ... rIdentNode  identifier for the node
void NuTo::StructureBase::GroupAddNode(const std::string& rIdentGroup, int rIdNode)
{
    boost::ptr_map<std::string,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNode] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=GroupBase::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNode] A node can be added only to a node group.");

    itGroup->second->AddMember(NodeGetNodePtr(rIdNode));
}

//! @brief ... Adds a node to a node group
//! @param ... rIdentGroup identifier for the group
//! @param ... rIdentNode  identifier for the node
void NuTo::StructureBase::GroupAddElement(const std::string& rIdentGroup, int rIdElement)
{
    boost::ptr_map<std::string,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddElement] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=GroupBase::Elements)
        throw MechanicsException("[NuTo::StructureBase::GroupAddElement] An element can be added only to an element group.");

    itGroup->second->AddMember(ElementGetElementPtr(rIdElement));
}

//! @brief ... Unites two groups and stores the result in a new group
//! @param ... rIdentGroup1 identifier for the first group
//! @param ... rIdentGroup2 identifier for the second group
//! @result ... rIdentGroupResult identifier for the created result group
void NuTo::StructureBase::GroupUnion(const std::string& rIdentGroup1, const std::string& rIdentGroup2, const std::string& rIdentGroupResult)
{
    boost::ptr_map<std::string,GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group1 with the given identifier does not exist.");
    boost::ptr_map<std::string,GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group2 with the given identifier does not exist.");
    boost::ptr_map<std::string,GroupBase>::iterator itGroupResult = mGroupMap.find(rIdentGroupResult);
    if (itGroupResult!=mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] GroupResult with the given identifier already exists.");

    // insert the member
    mGroupMap.insert(const_cast<std::string&>(rIdentGroupResult), (*itGroup1).second->Unite((*itGroup2).second));
}

//! @brief ... Difference between two groups and stores the result in a new group
//! @param ... rIdentGroup1 identifier for the first group
//! @param ... rIdentGroup2 identifier for the second group
//! @result ... rIdentGroupResult identifier for the created result group
void NuTo::StructureBase::GroupDifference(const std::string& rIdentGroup1, const std::string& rIdentGroup2, const std::string& rIdentGroupResult)
{
    boost::ptr_map<std::string,GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group1 with the given identifier does not exist.");
    boost::ptr_map<std::string,GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group2 with the given identifier does not exist.");
    boost::ptr_map<std::string,GroupBase>::iterator itGroupResult = mGroupMap.find(rIdentGroupResult);
    if (itGroupResult!=mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] GroupResult with the given identifier already exists.");

    // insert the member
    mGroupMap.insert(const_cast<std::string&>(rIdentGroupResult), (*itGroup1).second->Difference((*itGroup2).second));
}

//! @brief ... Calculates the intersection between two groups and stores the result in a new group
//! @param ... rIdentGroup1 identifier for the first group
//! @param ... rIdentGroup2 identifier for the second group
//! @result ... rIdentGroupResult identifier for the created result group
void NuTo::StructureBase::GroupIntersection(const std::string& rIdentGroup1, const std::string& rIdentGroup2, const std::string& rIdentGroupResult)
{
    boost::ptr_map<std::string,GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group1 with the given identifier does not exist.");
    boost::ptr_map<std::string,GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group2 with the given identifier does not exist.");
    boost::ptr_map<std::string,GroupBase>::iterator itGroupResult = mGroupMap.find(rIdentGroupResult);
    if (itGroupResult!=mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] GroupResult with the given identifier already exists.");

    // insert the member
    mGroupMap.insert(const_cast<std::string&>(rIdentGroupResult), (*itGroup1).second->Intersection((*itGroup2).second));
}

//! @brief ... Calculates the symmetric difference between two groups and stores the result in a new group
//! @param ... rIdentGroup1 identifier for the first group
//! @param ... rIdentGroup2 identifier for the second group
//! @result ... rIdentGroupResult identifier for the created result group
void NuTo::StructureBase::GroupSymmetricDifference(const std::string& rIdentGroup1, const std::string& rIdentGroup2, const std::string& rIdentGroupResult)
{
    boost::ptr_map<std::string,GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group1 with the given identifier does not exist.");
    boost::ptr_map<std::string,GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group2 with the given identifier does not exist.");
    boost::ptr_map<std::string,GroupBase>::iterator itGroupResult = mGroupMap.find(rIdentGroupResult);
    if (itGroupResult!=mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] GroupResult with the given identifier already exists.");

    // insert the member
    mGroupMap.insert(const_cast<std::string&>(rIdentGroupResult), (*itGroup1).second->SymmetricDifference((*itGroup2).second));
}

//! @brief ... Returns the number of members in a group
//! @param ... rIdentGroup identifier for the group
//! @return ... number of members
int NuTo::StructureBase::GroupGetNumMembers(const std::string& rIdentGroup)const
{
    boost::ptr_map<std::string,GroupBase>::const_iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupGetNumMembers] Group with the given identifier does not exist.");
    return (*itGroup).second->GetNumMembers();
}

// get group pointer from group identifier
NuTo::GroupBase* NuTo::StructureBase::GroupGetGroupPtr(const std::string& rIdent)
{
    boost::ptr_map<std::string,GroupBase>::iterator it = this->mGroupMap.find(rIdent);
    if (it == this->mGroupMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::GroupGetGroupPtr] Group does not exist.");
    }
    return it->second;
}

// get section pointer from section identifier
const NuTo::GroupBase* NuTo::StructureBase::GroupGetGroupPtr(const std::string& rIdent) const
{
    boost::ptr_map<std::string,GroupBase>::const_iterator it = this->mGroupMap.find(rIdent);
    if (it == this->mGroupMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::GroupGetGroupPtr] Group does not exist.");
    }
    return it->second;
}
