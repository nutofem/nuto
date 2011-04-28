// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"

//! @brief sets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeSetDisplacements(int rNode, const FullMatrix<double>& rDisplacements)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	NodeBase* nodePtr=NodeGetNodePtr(rNode);
	this->mUpdateTmpStaticDataRequired=true;

	if (rDisplacements.GetNumColumns()!=1)
	throw MechanicsException("[NuTo::StructureBase::NodeSetDisplacements] Displacement matrix has to have a single column.");
	try
	{
		switch (rDisplacements.GetNumRows())
		{
		case 1:
			nodePtr->SetDisplacements1D(rDisplacements.mEigenMatrix.data());
		break;
		case 2:
			nodePtr->SetDisplacements2D(rDisplacements.mEigenMatrix.data());
		break;
		case 3:
			nodePtr->SetDisplacements3D(rDisplacements.mEigenMatrix.data());
		break;
		default:
			throw MechanicsException("[NuTo::StructureBase::NodeSetDisplacements] The number of displacement components is either 1, 2 or 3.");
		}
	}
    catch(NuTo::MechanicsException & b)
	{
    	b.AddMessage("[NuTo::StructureBase::NodeSetDisplacements] Error setting displacements.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeSetDisplacements] Error setting displacements of node (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureBase::NodeSetDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief sets the displacements of a group of nodes
//! @param rIdent node group identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeGroupSetDisplacements(int rGroupIdent, const FullMatrix<double>& rDisplacements)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	this->mUpdateTmpStaticDataRequired=true;
	if (rDisplacements.GetNumColumns()!=1)
	     throw MechanicsException("[NuTo::StructureBase::NodeGroupSetDisplacements] Displacement matrix has to have a single column.");

	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::NodeGroupSetDisplacements] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
    	throw MechanicsException("[NuTo::StructureBase::NodeGroupSetDisplacements] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
		try
		{
			switch (rDisplacements.GetNumRows())
			{
			case 1:
				itNode->second->SetDisplacements1D(rDisplacements.mEigenMatrix.data());
			break;
			case 2:
				itNode->second->SetDisplacements2D(rDisplacements.mEigenMatrix.data());
			break;
			case 3:
				itNode->second->SetDisplacements3D(rDisplacements.mEigenMatrix.data());
			break;
			default:
				throw MechanicsException("[NuTo::StructureBase::NodeGroupSetDisplacements] The number of displacement components is either 1, 2 or 3.");
			}
		}
		catch(NuTo::MechanicsException & b)
		{
			b.AddMessage("[NuTo::StructureBase::NodeGroupSetDisplacements] Error setting displacements.");
			throw b;
		}
		catch(...)
		{
			throw MechanicsException("[NuTo::StructureBase::NodeGroupSetDisplacements] Error setting displacements of node (unspecified exception).");
		}
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::NodeGroupSetDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief gets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeGetDisplacements(int rNode, FullMatrix<double>& rDisplacements)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	const NodeBase* nodePtr = NodeGetNodePtr(rNode);

	try
	{
		switch (nodePtr->GetNumDisplacements())
		{
		case 1:
			rDisplacements.Resize(1,1);
			nodePtr->GetDisplacements1D(rDisplacements.mEigenMatrix.data());
		break;
		case 2:
			rDisplacements.Resize(2,1);
			nodePtr->GetDisplacements2D(rDisplacements.mEigenMatrix.data());
		break;
		case 3:
			rDisplacements.Resize(3,1);
			nodePtr->GetDisplacements3D(rDisplacements.mEigenMatrix.data());
		break;
		case 0:
			throw MechanicsException("[NuTo::StructureBase::NodeGetDisplacements] Node has no displacements.");
		break;
		}
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeGetDisplacements] Error getting displacements.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeGetDisplacements] Error getting displacements of node (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureBase::NodeGetDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief gets the coordinates of a node
//! @param rNode node identifier
//! @param rCoordinates matrix (one column) with the coordinates
void NuTo::StructureBase::NodeGetCoordinates(int rNode, NuTo::FullMatrix<double>& rCoordinates)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	const NodeBase* nodePtr = NodeGetNodePtr(rNode);

	try
	{
		switch (nodePtr->GetNumCoordinates())
		{
		case 1:
			rCoordinates.Resize(1,1);
			nodePtr->GetCoordinates1D(rCoordinates.mEigenMatrix.data());
		break;
		case 2:
			rCoordinates.Resize(2,1);
			nodePtr->GetCoordinates2D(rCoordinates.mEigenMatrix.data());
		break;
		case 3:
			rCoordinates.Resize(3,1);
			nodePtr->GetCoordinates3D(rCoordinates.mEigenMatrix.data());
		break;
		case 0:
			throw MechanicsException("[NuTo::StructureBase::NodeGetCoordinates] Node has no coordinates.");
		break;
		}
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeGetCoordinates] Error getting coordinates.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeGetCoordinates] Error getting coordinates of node (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureBase::NodeGetCoordinates] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief calculate the internal force vector for a node
//! @param rId ... node id
//! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
void NuTo::StructureBase::NodeInternalForce(int rId, NuTo::FullMatrix<double>& rNodeForce) const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	try
	{
        const NodeBase* nodePtr = NodeGetNodePtr(rId);
        NodeInternalForce(nodePtr,rNodeForce);
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeGradientInternalPotential] Error getting gradient of internal potential.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeGradientInternalPotential] Error getting gradient of internal potential (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::NodeInternalForce] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief calculate the internal force vector for a node group of nodes
//! @param rGroupIdent ... group identifier
//! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
void NuTo::StructureBase::NodeGroupInternalForce(int rGroupIdent, NuTo::FullMatrix<double>& rNodeForce) const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::NodeGroupForce] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
    	throw MechanicsException("[NuTo::StructureBase::NodeGroupForce] Group is not a node group.");
    const Group<NodeBase> *nodeGroup = dynamic_cast<const Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

	NuTo::FullMatrix<double>nodeForceLocal;

	if (nodeGroup->GetNumMembers()==0)
		throw MechanicsException("[NuTo::StructureBase::NodeGroupForce] Node group is empty.");
	rNodeForce.Resize(nodeGroup->begin()->second->GetNumDisplacements(),1);

    for (Group<NodeBase>::const_iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
		try
		{
			NodeInternalForce(itNode->second, nodeForceLocal);
			if (nodeForceLocal.GetNumRows()!=rNodeForce.GetNumRows())
				throw MechanicsException("[NuTo::StructureBase::NodeGroupForce] The number of displacement components is not equal for all members of the group.");
			rNodeForce+=nodeForceLocal;
		}
		catch(NuTo::MechanicsException & b)
		{
			b.AddMessage("[NuTo::StructureBase::NodeGroupForce] Error getting gradient of internal potential.");
			throw b;
		}
		catch(...)
		{
			throw MechanicsException("[NuTo::StructureBase::NodeGroupForce] Error getting gradient of internal potential (unspecified exception).");
		}
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::NodeGroupInternalForce] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief calculate the internal force vector for a node
//! @param rNodePtr  node for which this has to be calculated
//! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
void NuTo::StructureBase::NodeInternalForce(const NodeBase* rNodePtr, NuTo::FullMatrix<double>& rNodeForce) const
{
	try
	{
		rNodeForce.Resize(rNodePtr->GetNumDisplacements(),1);

		//go through all elements and check, if the node belongs to the element
		std::vector<const ElementBase*> elements;
		GetElementsTotal(elements);
		for (unsigned int countElement=0; countElement<elements.size(); countElement++)
		{
			const ElementBase* elementPtr=elements[countElement];
			for (int countNode=0; countNode<elementPtr->GetNumNodes(); countNode++)
			{
				if (elementPtr->GetNode(countNode)==rNodePtr)
				{
					NuTo::FullMatrix<double> result;
					std::vector<int> globalDofs;
					elementPtr->CalculateGradientInternalPotential(result,globalDofs);

					for (int countDof=0; countDof< rNodePtr->GetNumDisplacements(); countDof++)
					{
                        int theDof = rNodePtr->GetDofDisplacement(countDof);
                        for (unsigned int countGlobalDofs=0; countGlobalDofs<globalDofs.size(); countGlobalDofs++)
                        {
                        	if (globalDofs[countGlobalDofs] == theDof)
                        	{
                        		rNodeForce(countDof,0)+=result(countGlobalDofs,0);
                        	}
                        }
					}
				}
			}
		}
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeGradientInternalPotential] Error getting gradient of internal potential.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeGradientInternalPotential] Error getting gradient of internal potential (unspecified exception).");
	}

}

//! @brief ... store all element ids connected to this node in a vector
//! @param rNode (Input) 			... node id
//! @param rElementNumbers (Output) ... vector of element ids
void NuTo::StructureBase::NodeGetElements(const int rNodeId, NuTo::FullMatrix<int>& rElementNumbers)
{
    throw MechanicsException("[NuTo::StructureBase::NodeGetElements] Not available for this structure type.");
}

//! @brief ... store all elements connected to this node in a vector
//! @param rNode (Input) 		... node pointer
//! @param rElements (Output) 	... vector of element pointers
void NuTo::StructureBase::NodeGetElements(const NuTo::NodeBase* rNodePtr, std::vector<NuTo::ElementBase*>& rElements)
{
    throw MechanicsException("[NuTo::StructureBase::NodeGetElements] Not available for this structure type.");
}
