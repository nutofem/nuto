
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
BOOST_CLASS_EXPORT_GUID(NuTo::Group<NuTo::NodeBase>, "NuToGroupNodeBase")
BOOST_CLASS_EXPORT_GUID(std::set<NuTo::NodeBase*>, "stdSetNodeBase")
BOOST_CLASS_EXPORT_GUID(NuTo::Group<NuTo::ElementBase>, "NuToGroupElementBase")
BOOST_CLASS_EXPORT_GUID(std::set<NuTo::ElementBase*>, "stdSetElementBase")
#endif // ENABLE_SERIALIZATION
