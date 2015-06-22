#include "nuto/mechanics/structures/unstructured/Structure.h"

//! @brief ... Adds all elements to a group based on the type
//! @param ... rIdentGroup identifier for the group
//! @param ... rInterpolationType  identifier for the interpolation type
void NuTo::Structure::GroupAddElementFromType(int rIdentGroup, int rInterpolationType)
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
			throw MechanicsException("[NuTo::Structure::GroupAddElementFromType] Group is not an element group.");

		boost::ptr_map<int,InterpolationType>::iterator itInterpolationType = mInterpolationTypeMap.find(rInterpolationType);
		if (itInterpolationType==mInterpolationTypeMap.end())
            throw MechanicsException("[NuTo::Structure::GroupAddElementFromType] InterpolationType with the given identifier does not exist.");

		InterpolationType* interpolationType = itInterpolationType->second;

		for (boost::ptr_map<int,ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
	    {
		    ElementBase* element = elementIter->second;

		    if (element->GetInterpolationType() == interpolationType)
				itGroup->second->AddMember(elementIter->first, elementIter->second);
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
