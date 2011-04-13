// $Id$

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"


namespace NuTo
{
template<>
Groups::eGroupId Group<NodeBase>::GetType()const
{
    return Groups::Nodes;
}

template<>
Groups::eGroupId Group<ElementBase>::GetType()const
{
    return Groups::Elements;
}
template<>
std::string Group<NodeBase>::GetTypeString()const
{
    return std::string("Nodes");
}

template<>
std::string Group<ElementBase>::GetTypeString()const
{
    return std::string("Elements");
}

//! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
template<>
Group<ElementBase>* Group<ElementBase>::AsGroupElement()
{
    return this;
}

//! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
template<>
const Group<ElementBase>* Group<ElementBase>::AsGroupElement()const
{
    return this;
}

//! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
template<>
Group<NodeBase>* Group<ElementBase>::AsGroupNode()
{
    throw MechanicsException("[Group<ElementBase>::AsGroupNode] group is not a node group");
}

//! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
template<>
const Group<NodeBase>* Group<ElementBase>::AsGroupNode()const
{
    throw MechanicsException("[Group<ElementBase>::AsGroupNode] group is not a node group");
}

//! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
template<>
Group<ElementBase>* Group<NodeBase>::AsGroupElement()
{
    throw MechanicsException("[Group<NodeBase>::AsGroupNode] group is not an element group");
}

//! @brief either casts the pointer to an element group or throws an exception for groups which are not element groups
template<>
const Group<ElementBase>* Group<NodeBase>::AsGroupElement()const
{
    throw MechanicsException("[Group<NodeBase>::AsGroupNode] group is not an element group");
}

//! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
template<>
Group<NodeBase>* Group<NodeBase>::AsGroupNode()
{
    return this;
}

//! @brief either casts the pointer to a node group or throws an exception for groups which are not node groups
template<>
const Group<NodeBase>* Group<NodeBase>::AsGroupNode()const
{
    return this;
}

//! @brief info for group members
//! @param rVerboseLevel verbose Level
//! @param rStructure Structure that holds the members
template<>
void Group<NodeBase>::Info(int rVerboseLevel, const NuTo::StructureBase* rStructure)const
{
	std::cout << "    Number of members : " << this->size() << std::endl;
    if (rVerboseLevel>2)
    {
    	std::cout << "    members :" <<std::endl;
    	for (Group<NodeBase>::const_iterator it=this->begin(); it!= this->end(); it++)
    	{
    		std::cout << "      " << rStructure->NodeGetId(*it) << std::endl;
    	}
    }
}

//! @brief info for group members
//! @param rVerboseLevel verbose Level
//! @param rStructure Structure that holds the members
template<>
void Group<ElementBase>::Info(int rVerboseLevel, const NuTo::StructureBase* rStructure)const
{
	std::cout << "    Number of members : " << this->size() << std::endl;
    if (rVerboseLevel>2)
    {
    	std::cout << "    members :" <<std::endl;
    	for (Group<ElementBase>::const_iterator it=this->begin(); it!= this->end(); it++)
    	{
    		std::cout << "      " << rStructure->ElementGetId(*it) << std::endl;
    	}
    }
}

}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Group<NuTo::NodeBase>)
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Group<NuTo::ElementBase>)
#endif // ENABLE_SERIALIZATION
