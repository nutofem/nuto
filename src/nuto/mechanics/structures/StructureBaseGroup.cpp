// $Id$

#include <set>

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"

void NuTo::StructureBase::GroupInfo(int rVerboseLevel)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    std::cout << "number of groups  : " << mGroupMap.size() << std::endl;
    if (rVerboseLevel>2)
    {
        for (boost::ptr_map<int,GroupBase>::const_iterator it= mGroupMap.begin(); it!=mGroupMap.end(); it++)
        {
            std::cout << "  Group " << it->first << std::endl;
            it->second->Info(rVerboseLevel, this);
            std::cout <<  std::endl;
        }
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupInfo] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

}

//! @brief gives the identifier of a group
//! @param pointer to a group
//! @return identifier
int NuTo::StructureBase::GroupGetId(GroupBase* rGroup)const
{
    for (boost::ptr_map<int,GroupBase>::const_iterator
            it = mGroupMap.begin(); it!= mGroupMap.end(); it++)
    {
        if (it->second==rGroup)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::GroupGetId] Group does not exist.");
}

//! @brief ... Creates a group for the structure
//! @param ... rType  type of the group, e.g. "NODES" or "ELEMENTS"
//! @return ... rIdent identifier for the group
int NuTo::StructureBase::GroupCreate(const std::string& rType)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	//find unused integer id
	int groupNumber(mGroupMap.size());
	boost::ptr_map<int,GroupBase>::iterator it = mGroupMap.find(groupNumber);
	while (it!=mGroupMap.end())
	{
		groupNumber++;
		it = mGroupMap.find(groupNumber);
	}

	// transform string to uppercase
    std::string GroupTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(GroupTypeString), (int(*)(int)) toupper);

    if (GroupTypeString==std::string("ELEMENTS"))
    {
        GroupCreate(groupNumber,NuTo::Groups::Elements);
    }
    else
    {
		if (GroupTypeString==std::string("NODES"))
		{
			GroupCreate(groupNumber,NuTo::Groups::Nodes);
		}
		else
		{
			throw MechanicsException("[NuTo::StructureBase::GroupCreate] Group type not implemented.");
		}
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupCreate] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

    return groupNumber;
}

//! @brief ... Creates a group for the structure
//! @param ... rIdent identifier for the group
//! @param ... rType  type of the group
void NuTo::StructureBase::GroupCreate(int id, NuTo::Groups::eGroupId rEnumType)
{
	boost::ptr_map<int,GroupBase>::iterator it = mGroupMap.find(id);
	if (it!=mGroupMap.end())
		throw MechanicsException("[NuTo::StructureBase::GroupCreate] Group id already exists.");

    switch(rEnumType)
    {
    case NuTo::Groups::Elements:
    	mGroupMap.insert(id, new NuTo::Group<ElementBase>);
    	break;
    case NuTo::Groups::Nodes:
    	mGroupMap.insert(id, new NuTo::Group<NodeBase>);
    	break;
    default:
    	throw MechanicsException("[NuTo::StructureBase::GroupCreate] Group type not implemented.");
    }
}

//! @brief ... Deletes a group from the structure
//! @param ... rIdent identifier for the group
void NuTo::StructureBase::GroupDelete(int rIdentGroup)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // find group in map
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup == this->mGroupMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::GroupDelete] Group with the given identifier does not exist.");
    }
    else
    {
        // delete load from map
        this->mGroupMap.erase(itGroup);
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupDelete] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief ... Adds a node to a node group
//! @param ... rIdentGroup identifier for the group
//! @param ... rIdentNode  identifier for the node
void NuTo::StructureBase::GroupAddNode(int rIdentGroup, int rIdNode)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNode] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNode] A node can be added only to a node group.");

    itGroup->second->AddMember(rIdNode, NodeGetNodePtr(rIdNode));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupAddNode] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief ... Adds an element to an element group
//! @param ... rIdentGroup identifier for the group
//! @param ... rIdentNode  identifier for the element
void NuTo::StructureBase::GroupAddElement(int rIdentGroup, int rIdElement)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddElement] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::GroupAddElement] An element can be added only to an element group.");

    itGroup->second->AddMember(rIdElement, ElementGetElementPtr(rIdElement));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupAddElement] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief ... Adds all nodes to a group whose coordinates are in the specified range
//! @param ... rIdentGroup identifier for the group
//! @param ... rDirection either 0,1,2 for x,y, or z
//! @param ... rMin ... minimum value
//! @param ... rMax ... maximum value
void NuTo::StructureBase::GroupAddNodeCoordinateRange(int rIdentGroup, int rDirection, double rMin, double rMax)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
   boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] A node can be added only to a node group.");

    if(rDirection<0 || rDirection>mDimension)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] The direction is either 0(x),1(Y) or 2(Z) and has to be smaller than the dimension of the structure.");

    std::vector<std::pair<int,NodeBase*> > nodeVector;
    this->GetNodesTotal(nodeVector);
    double coordinates[3];

    for (unsigned int countNode=0; countNode<nodeVector.size(); countNode++)
    {
    	NodeBase* nodePtr(nodeVector[countNode].second);
    	if (nodePtr->GetNumCoordinates()<1)
    		continue;
    	switch (mDimension)
    	{
    	case 1:
    	    nodePtr->GetCoordinates1D(coordinates);
    	    break;
    	case 2:
    	    nodePtr->GetCoordinates2D(coordinates);
    	    break;
    	case 3:
    	    nodePtr->GetCoordinates3D(coordinates);
    	    break;
    	default:
    		throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCoordinateRange] unsupported dimension of the structure.");
    	}
    	if (coordinates[rDirection]>=rMin && coordinates[rDirection]<=rMax)
            itGroup->second->AddMember(nodeVector[countNode].first,nodePtr);
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupAddNodeCoordinateRange] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief ... Adds all nodes to a group whose coordinates are in the specified range
//! @param ... rIdentGroup identifier for the group
//! @param ... rCenter center of the selection circle
//! @param ... rMin ... minimum radius
//! @param ... rMax ... maximum radius
void NuTo::StructureBase::GroupAddNodeRadiusRange(int rIdentGroup, NuTo::FullVector<double,Eigen::Dynamic> rCenter, double rMin, double rMax)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
   boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] A node can be added only to a node group.");

    if (rCenter.GetNumRows()!=mDimension || rCenter.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] The center point must have the same number of coordinates as the dimension of the structure.");

    if(rMin>rMax)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] The minimum radius must not be larger than the maximum radius.");

    std::vector<std::pair<int,NodeBase*> > nodeVector;
    this->GetNodesTotal(nodeVector);
    double coordinates[3];
    double rMin2 = rMin*rMin;
    double rMax2 = rMax*rMax;

    for (unsigned int countNode=0; countNode<nodeVector.size(); countNode++)
    {
        NodeBase* nodePtr(nodeVector[countNode].second);
        if (nodePtr->GetNumCoordinates()<1)
            continue;
        double r2(0.);
        switch (mDimension)
        {
        case 1:
            nodePtr->GetCoordinates1D(coordinates);
            r2 = (coordinates[0]-rCenter(0))*(coordinates[0]-rCenter(0));
            break;
        case 2:
            nodePtr->GetCoordinates2D(coordinates);
            r2 = (coordinates[0]-rCenter(0))*(coordinates[0]-rCenter(0))
                +(coordinates[1]-rCenter(1))*(coordinates[1]-rCenter(1));
            break;
        case 3:
            nodePtr->GetCoordinates3D(coordinates);
            r2 = (coordinates[0]-rCenter(0))*(coordinates[0]-rCenter(0))
                +(coordinates[1]-rCenter(1))*(coordinates[1]-rCenter(1))
                +(coordinates[2]-rCenter(2))*(coordinates[2]-rCenter(2));
            break;
        default:
            throw MechanicsException("[NuTo::StructureBase::GroupAddNodeRadiusRange] unsupported dimension of the structure.");
        }
        if (r2>=rMin2 && r2<=rMax2)
            itGroup->second->AddMember(nodeVector[countNode].first,nodePtr);
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupAddNodeRadiusRange] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief ... Adds all nodes to a group whose coordinates are on a cylinder with the radius in the in the specified range
//! @param ... rIdentGroup identifier for the group
//! @param ... rCenter center of the cylinder
//! @param ... rAxis axis of the cylinder
//! @param ... rMin ... minimum radius
//! @param ... rMax ... maximum radius
void NuTo::StructureBase::GroupAddNodeCylinderRadiusRange(int rIdentGroup, NuTo::FullVector<double,Eigen::Dynamic> rCenter,
		 NuTo::FullVector<double,Eigen::Dynamic> rDirection, double rMin, double rMax)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
   boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Nodes)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] A node can be added only to a node group.");

    if (rCenter.GetNumRows()!=mDimension || rCenter.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] The center point must have the same number of coordinates as the dimension of the structure.");

    if (rDirection.GetNumRows()!=mDimension || rDirection.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] The direction point must have the same number of coordinates as the dimension of the structure.");

    if(rMin>rMax)
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] The minimum radius must not be larger than the maximum radius.");

    switch (mDimension)
    {
    case 2:
    {
        std::vector<std::pair<int,NodeBase*> > nodeVector;
        this->GetNodesTotal(nodeVector);
        FullVector<double,2> coordinates;
        FullVector<double,2> vecDelta;
        double rMin2 = rMin*rMin;
        double rMax2 = rMax*rMax;

        for (unsigned int countNode=0; countNode<nodeVector.size(); countNode++)
        {
            NodeBase* nodePtr(nodeVector[countNode].second);
            if (nodePtr->GetNumCoordinates()!=2)
                continue;
            double r2(0.);
    		nodePtr->GetCoordinates2D(coordinates.data());
    		vecDelta = coordinates - rCenter;

    		r2 = (vecDelta(0)*vecDelta(0) + vecDelta(1)*vecDelta(1));

    		if (r2>=rMin2 && r2<=rMax2)
    		{
                itGroup->second->AddMember(nodeVector[countNode].first,nodePtr);
    		}
        }
    	break;
    }
    case 3:
    {
        std::vector<std::pair<int,NodeBase*> > nodeVector;
        this->GetNodesTotal(nodeVector);
        FullVector<double,3> coordinates;
        FullVector<double,3> vecPtrCenter;
        FullVector<double,3> vecPtrProjection;
        FullVector<double,3> vecDelta;
        double rMin2 = rMin*rMin;
        double rMax2 = rMax*rMax;

        //normalize Diretion Vector
        rDirection*=1./rDirection.Norm();

        for (unsigned int countNode=0; countNode<nodeVector.size(); countNode++)
        {
            NodeBase* nodePtr(nodeVector[countNode].second);
            if (nodePtr->GetNumCoordinates()!=3)
                continue;
            double r2(0.);
    		nodePtr->GetCoordinates3D(coordinates.data());
    		vecPtrCenter = coordinates - rCenter;

    		//get projection onto axis
    		double s = rDirection.transpose()*vecPtrCenter;
    		vecPtrProjection = rCenter+rDirection*s;
    		vecDelta = coordinates - vecPtrProjection;

    		r2 = (vecDelta(0)*vecDelta(0) + vecDelta(1)*vecDelta(1) + vecDelta(2)*vecDelta(2));

    		if (r2>=rMin2 && r2<=rMax2)
    		{
                itGroup->second->AddMember(nodeVector[countNode].first,nodePtr);
    		}
        }
    	break;
    }
    default:
        throw MechanicsException("[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] Only implemented for 2D and 3D.");

    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupAddNodeCylinderRadiusRange] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief ... Adds all elements to a group whose nodes are in the given node group
//! @param ... rElementGroupId identifier for the element group
//! @param ... rNodeGroupId idenbtifier for the node group
//! @param ... rHaveAllNodes if set to true, the element is only selected when all element nodes are in the node group, if set
//! to false, the element is select if at least one node is in the node group
void NuTo::StructureBase::GroupAddElementsFromNodes(int rElementGroupId, int rNodeGroupId, bool rHaveAllNodes)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    //get element group
	Group<ElementBase> *elementGroup = this->GroupGetGroupPtr(rElementGroupId)->AsGroupElement();

    //get node group
	const Group<NodeBase> *nodeGroup = this->GroupGetGroupPtr(rNodeGroupId)->AsGroupNode();

	//since the search is done via the id's, the element nodes are ptr, so make another set with the node ptrs
	std::set<const NodeBase*> nodePtrSet;
    for (Group<NodeBase>::const_iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
    	nodePtrSet.insert(itNode->second);
    }

    std::vector<std::pair<int,ElementBase*> > elementVector;
    this->GetElementsTotal(elementVector);
	std::vector<const NodeBase*> elementNodes;
    for (unsigned int countElement=0; countElement<elementVector.size(); countElement++)
    {
        try
        {
        	if (!elementGroup->Contain(elementVector[countElement].first))
        	{
        		bool addElement;
        		if (rHaveAllNodes==true)
        		{
        			addElement = true;
					int countNode;
					for (countNode=0; (countNode<elementVector[countElement].second->GetNumNodes()) && (addElement==true);
							countNode++)
					{
						if (nodePtrSet.find(elementVector[countElement].second->GetNode(countNode))==nodePtrSet.end())
						{
							addElement=false;
						}
					}
        		}
        		else
        		{
        			addElement = false;
					int countNode;
					for (countNode=0; (countNode<elementVector[countElement].second->GetNumNodes()) && (addElement==false);
							countNode++)
					{
						if (nodePtrSet.find(elementVector[countElement].second->GetNode(countNode))!=nodePtrSet.end())
						{
							addElement=true;
						}
					}
        		}
				if (addElement)
				{
					//add the element;
					elementGroup->AddMember(elementVector[countElement].first, elementVector[countElement].second);
				}
        	}
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(elementVector[countElement].second)==elementVector[countElement].first);
            ss << elementVector[countElement].first;
            e.AddMessage("[NuTo::StructureBase::GroupAddElementsFromNodes] Error for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(elementVector[countElement].second)==elementVector[countElement].first);
            ss << elementVector[countElement].first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::GroupAddElementsFromNodes] Error for element " + ss.str() + ".");
        }
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupAddElementsFromNodes] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief ... Unites two groups and stores the result in a new group
//! @param ... rIdentGroup1 identifier for the first group
//! @param ... rIdentGroup2 identifier for the second group
//! @result ... rIdentGroupResult identifier for the created result group
int NuTo::StructureBase::GroupUnion(int rIdentGroup1, int rIdentGroup2)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group1 with the given identifier does not exist.");
    boost::ptr_map<int,GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupUnite] Group2 with the given identifier does not exist.");

    //find unused integer id
	int groupNumber(mGroupMap.size());
	boost::ptr_map<int,GroupBase>::iterator it = mGroupMap.find(groupNumber);
	while (it!=mGroupMap.end())
	{
		groupNumber++;
		it = mGroupMap.find(groupNumber);
	}
    // insert the member
    mGroupMap.insert(groupNumber, (*itGroup1).second->Unite((*itGroup2).second));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupUnite] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
    return groupNumber;
}

//! @brief ... Difference between two groups and stores the result in a new group
//! @param ... rIdentGroup1 identifier for the first group
//! @param ... rIdentGroup2 identifier for the second group
//! @result ... rIdentGroupResult identifier for the created result group
int NuTo::StructureBase::GroupDifference(int rIdentGroup1, int rIdentGroup2)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupDifference] Group1 with the given identifier does not exist.");
    boost::ptr_map<int,GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupDifference] Group2 with the given identifier does not exist.");
    //find unused integer id
	int groupNumber(mGroupMap.size());
	boost::ptr_map<int,GroupBase>::iterator it = mGroupMap.find(groupNumber);
	while (it!=mGroupMap.end())
	{
		groupNumber++;
		it = mGroupMap.find(groupNumber);
	}

    // insert the member
    mGroupMap.insert(groupNumber, (*itGroup1).second->Difference((*itGroup2).second));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupDifference] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

    return groupNumber;
}

//! @brief ... Calculates the intersection between two groups and stores the result in a new group
//! @param ... rIdentGroup1 identifier for the first group
//! @param ... rIdentGroup2 identifier for the second group
//! @result ... rIdentGroupResult identifier for the created result group
int NuTo::StructureBase::GroupIntersection(int rIdentGroup1, int rIdentGroup2)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupIntersection] Group1 with the given identifier does not exist.");
    boost::ptr_map<int,GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupIntersection] Group2 with the given identifier does not exist.");

    //find unused integer id
	int groupNumber(mGroupMap.size());
	boost::ptr_map<int,GroupBase>::iterator it = mGroupMap.find(groupNumber);
	while (it!=mGroupMap.end())
	{
		groupNumber++;
		it = mGroupMap.find(groupNumber);
	}

    // insert the member
    mGroupMap.insert(groupNumber, (*itGroup1).second->Intersection((*itGroup2).second));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupIntersection] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

    return groupNumber;
}

//! @brief ... Calculates the symmetric difference between two groups and stores the result in a new group
//! @param ... rIdentGroup1 identifier for the first group
//! @param ... rIdentGroup2 identifier for the second group
//! @result ... rIdentGroupResult identifier for the created result group
int NuTo::StructureBase::GroupSymmetricDifference(int rIdentGroup1, int rIdentGroup2)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup1 = mGroupMap.find(rIdentGroup1);
    if (itGroup1==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupSymmetricDifference] Group1 with the given identifier does not exist.");
    boost::ptr_map<int,GroupBase>::iterator itGroup2 = mGroupMap.find(rIdentGroup2);
    if (itGroup2==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupSymmetricDifference] Group2 with the given identifier does not exist.");

    //find unused integer id
	int groupNumber(mGroupMap.size());
	boost::ptr_map<int,GroupBase>::iterator it = mGroupMap.find(groupNumber);
	while (it!=mGroupMap.end())
	{
		groupNumber++;
		it = mGroupMap.find(groupNumber);
	}

    // insert the member
    mGroupMap.insert(groupNumber, (*itGroup1).second->SymmetricDifference((*itGroup2).second));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupSymmetricDifference] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

    return groupNumber;
}

//! @brief ... Returns the number of members in a group
//! @param ... rIdentGroup identifier for the group
//! @return ... number of members
int NuTo::StructureBase::GroupGetNumMembers(int rIdentGroup)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupGetNumMembers] Group with the given identifier does not exist.");
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupGetNumMembers] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
    return (*itGroup).second->GetNumMembers();
}

//! @brief ... Returns a vector with the members of a group
//! @param ... rIdentGroup identifier for the group
//! @return ... vector of members
NuTo::FullVector<int, Eigen::Dynamic> NuTo::StructureBase::GroupGetMemberIds(int rIdentGroup)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupGetNumMembers] Group with the given identifier does not exist.");
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupGetNumMembers] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
    return (*itGroup).second->GetMemberIds();
}


//! @brief ... checks for a member in a group
//! @param ... rIdentGroup identifier for the group
//! @return ... rMember id (element id, node id etc.)
bool NuTo::StructureBase::GroupContainsMember(int rIdentGroup, int rMember)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::GroupContains] Group with the given identifier does not exist.");
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::GroupContains] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
    return (*itGroup).second->Contain(rMember);
}

// get group pointer from group identifier
NuTo::GroupBase* NuTo::StructureBase::GroupGetGroupPtr(int rIdent)
{
    boost::ptr_map<int,GroupBase>::iterator it = this->mGroupMap.find(rIdent);
    if (it == this->mGroupMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::GroupGetGroupPtr] Group does not exist.");
    }
    return it->second;
}

// get section pointer from section identifier
const NuTo::GroupBase* NuTo::StructureBase::GroupGetGroupPtr(int rIdent) const
{
    boost::ptr_map<int,GroupBase>::const_iterator it = this->mGroupMap.find(rIdent);
    if (it == this->mGroupMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::GroupGetGroupPtr] Group does not exist.");
    }
    return it->second;
}
