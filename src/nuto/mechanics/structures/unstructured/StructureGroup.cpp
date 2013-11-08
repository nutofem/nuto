#include "nuto/mechanics/structures/unstructured/Structure.h"

//! @brief ... Adds all elements to a group based on the type
//! @param ... rIdentGroup identifier for the group
//! @param ... rElemTypeStr  identifier for the element type
void NuTo::Structure::GroupAddElementFromType(int rIdentGroup, std::string rElemTypeStr)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
		boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
		if (itGroup==mGroupMap.end())
			throw MechanicsException("[NuTo::Structure::GroupAddElementFromType] Group with the given identifier does not exist.");
		if (itGroup->second->GetType()!=Groups::Elements)
			throw MechanicsException("[NuTo::Structure::GroupAddElementFromType] An element can be added only to an element group.");

		NuTo::Element::eElementType elementType = this->ElementTypeGetEnum(rElemTypeStr);

		boost::ptr_map<int,ElementBase>::iterator ElementIter = this->mElementMap.begin();
	    while (ElementIter != this->mElementMap.end())
	    {
			if (ElementIter->second->GetEnumType()==elementType)
				itGroup->second->AddMember(ElementIter->first, ElementIter->second);
	    }
    }
    catch(NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::Structure::GroupAddElementFromType] Error adding element.");
        throw e;
    }
    catch(...)
    {
        throw NuTo::MechanicsException
           ("[NuTo::Structure::GroupAddElementFromType] Error adding element.");
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::GroupAddElementFromType] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief ... Adds a node to a node group
//! @param ... rIdentGroup identifier for the group
//! @param ... rIdentNode  identifier for the node
void NuTo::Structure::GroupAddElement(int rIdentGroup, int rIdElement)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::GroupAddElement] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Elements)
        throw MechanicsException("[NuTo::Structure::GroupAddElement] An element can be added only to an element group.");

    itGroup->second->AddMember(rIdElement, ElementGetElementPtr(rIdElement));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::GroupAddElement] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}
