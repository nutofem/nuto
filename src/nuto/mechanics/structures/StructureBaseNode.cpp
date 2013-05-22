// $Id$

#include <boost/assign/ptr_map_inserter.hpp>

#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputFullVectorDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"

//! @brief sets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeSetDisplacements(int rNode, const FullVector<double,Eigen::Dynamic>& rDisplacements)
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
			nodePtr->SetDisplacements1D(rDisplacements.data());
		break;
		case 2:
			nodePtr->SetDisplacements2D(rDisplacements.data());
		break;
		case 3:
			nodePtr->SetDisplacements3D(rDisplacements.data());
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

//! @brief sets the displacements of a node
//! @param rIdent node identifier
//! @param rTimeDerivative time derivative (0 disp, 1 vel, 2 acc)
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeSetDisplacements(int rNode, int rTimeDerivative, const FullVector<double,Eigen::Dynamic>& rDisplacements)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	NodeBase* nodePtr=NodeGetNodePtr(rNode);
	this->mUpdateTmpStaticDataRequired=true;

	if (rDisplacements.GetNumColumns()!=1)
	    throw MechanicsException("[NuTo::StructureBase::NodeSetDisplacements] Displacement matrix has to have a single column.");
	if (nodePtr->GetNumTimeDerivatives()>rTimeDerivative)
	    throw MechanicsException("[NuTo::StructureBase::NodeSetDisplacements] number of time derivatives stored at node is less than the required value.");
	try
	{
		switch (rDisplacements.GetNumRows())
		{
		case 1:
			nodePtr->SetDisplacements1D(rTimeDerivative,rDisplacements.data());
		break;
		case 2:
			nodePtr->SetDisplacements2D(rTimeDerivative,rDisplacements.data());
		break;
		case 3:
			nodePtr->SetDisplacements3D(rTimeDerivative,rDisplacements.data());
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

//! @brief sets the rotations of a node
//! @param rIdent node identifier
//! @param rRotations matrix (one column) with the rotations
void NuTo::StructureBase::NodeSetRotations(int rNode, const FullVector<double,Eigen::Dynamic>& rRotations)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	NodeBase* nodePtr=NodeGetNodePtr(rNode);
	this->mUpdateTmpStaticDataRequired=true;

	if (rRotations.GetNumColumns()!=1)
	throw MechanicsException("[NuTo::StructureBase::NodeSetRotations] rotation matrix has to have a single column.");
	try
	{
		switch (rRotations.GetNumRows())
		{
		case 1:
			nodePtr->SetRotations2D(rRotations.data());
		break;
		case 3:
			nodePtr->SetRotations3D(rRotations.data());
		break;
		default:
			throw MechanicsException("[NuTo::StructureBase::NodeSetRotations] The number of rotation components is either 1, 3.");
		}
	}
    catch(NuTo::MechanicsException & b)
	{
    	b.AddMessage("[NuTo::StructureBase::NodeSetRotations] Error setting rotations.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeSetRotations] Error setting rotations of node (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureBase::NodeSetRotations] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief sets the displacements of a group of nodes
//! @param rIdent node group identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeGroupSetDisplacements(int rGroupIdent, const FullVector<double,Eigen::Dynamic>& rDisplacements)
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
				itNode->second->SetDisplacements1D(rDisplacements.data());
			break;
			case 2:
				itNode->second->SetDisplacements2D(rDisplacements.data());
			break;
			case 3:
				itNode->second->SetDisplacements3D(rDisplacements.data());
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

//! @brief sets the displacements of a group of nodes
//! @param rIdent node group identifier
//! @param rTimeDerivative time derivative (0 disp, 1 vel, 2 acc)
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureBase::NodeGroupSetDisplacements(int rGroupIdent, int rTimeDerivative, const FullVector<double,Eigen::Dynamic>& rDisplacements)
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
			if (itNode->second->GetNumTimeDerivatives()<rTimeDerivative)
                throw MechanicsException("[NuTo::StructureBase::NodeGroupSetDisplacements] not does not have a sufficient number of time derivatives.");

			switch (rDisplacements.GetNumRows())
			{
			case 1:
				itNode->second->SetDisplacements1D(rTimeDerivative, rDisplacements.data());
			break;
			case 2:
				itNode->second->SetDisplacements2D(rTimeDerivative, rDisplacements.data());
			break;
			case 3:
				itNode->second->SetDisplacements3D(rTimeDerivative, rDisplacements.data());
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
void NuTo::StructureBase::NodeGetDisplacements(int rNode, FullVector<double,Eigen::Dynamic>& rDisplacements)const
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
			rDisplacements.Resize(1);
			nodePtr->GetDisplacements1D(rDisplacements.data());
		break;
		case 2:
			rDisplacements.Resize(2);
			nodePtr->GetDisplacements2D(rDisplacements.data());
		break;
		case 3:
			rDisplacements.Resize(3);
			nodePtr->GetDisplacements3D(rDisplacements.data());
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

//! @brief gets the rotations of a node
//! @param rIdent node identifier
//! @param rRotation matrix (one column) with the rotations
void NuTo::StructureBase::NodeGetRotations(int rNode, FullVector<double,Eigen::Dynamic>& rRotations)const
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
		case 2:
			rRotations.Resize(1);
			nodePtr->GetRotations2D(rRotations.data());
		break;
		case 3:
			rRotations.Resize(3);
			nodePtr->GetRotations3D(rRotations.data());
		break;
		default:
			throw MechanicsException("[NuTo::StructureBase::NodeGetRotations] Node has neither 1(2D) or 3(3D) rotations.");
		break;
		}
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeGetRotations] Error getting rotations.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeGetRotations] Error getting rotations of node (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureBase::NodeGetRotations] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}
//! @brief gets the displacements of a group of nodes
//! @param rNodeGroup node group identifier
//! @param rDisplacements matrix (rows/nodes columns/rDisplacements)
void NuTo::StructureBase::NodeGroupGetDisplacements(int rGroupIdent, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDisplacements)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::NodeGroupGetDisplacements] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
    	throw MechanicsException("[NuTo::StructureBase::NodeGroupGetDisplacements] Group is not a node group.");
    Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();
    assert(nodeGroup!=0);

    //all nodes have to have the same dimension
    if(nodeGroup->GetNumMembers()<1)
    	throw MechanicsException("[NuTo::StructureBase::NodeGroupGetDisplacements] Group has no members.");

    int numDisp= nodeGroup->begin()->second->GetNumDisplacements();
    //resize the matrix
    rDisplacements.Resize(nodeGroup->GetNumMembers(),numDisp);
	double disp[3];
	int theNode(0);
    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++, theNode++)
    {
		try
		{
			switch (numDisp)
			{
			case 1:
				itNode->second->GetDisplacements1D(disp);
				rDisplacements(theNode,0)=disp[0];
			break;
			case 2:
				itNode->second->GetDisplacements2D(disp);
				rDisplacements(theNode,0)=disp[0];
				rDisplacements(theNode,1)=disp[1];
			break;
			case 3:
				itNode->second->GetDisplacements3D(disp);
				rDisplacements(theNode,0)=disp[0];
				rDisplacements(theNode,1)=disp[1];
				rDisplacements(theNode,2)=disp[2];
			break;
			default:
				throw MechanicsException("[NuTo::StructureBase::NodeGroupGetDisplacements] The number of displacement components is either 1, 2 or 3.");
			}
		}
		catch(NuTo::MechanicsException & b)
		{
			b.AddMessage("[NuTo::StructureBase::NodeGroupGetDisplacements] Error getting displacements.");
			throw b;
		}
		catch(...)
		{
			throw MechanicsException("[NuTo::StructureBase::NodeGroupGetDisplacements] Error getting displacements of node (unspecified exception).");
		}
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::NodeGroupGetDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief gets the coordinates of a node
//! @param rNode node identifier
//! @param rCoordinates matrix (one column) with the coordinates
void NuTo::StructureBase::NodeGetCoordinates(int rNode, NuTo::FullVector<double,Eigen::Dynamic>& rCoordinates)const
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
			rCoordinates.Resize(1);
			nodePtr->GetCoordinates1D(rCoordinates.data());
		break;
		case 2:
			rCoordinates.Resize(2);
			nodePtr->GetCoordinates2D(rCoordinates.data());
		break;
		case 3:
			rCoordinates.Resize(3);
			nodePtr->GetCoordinates3D(rCoordinates.data());
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

//! @brief gets the coordinates of a group of nodes
//! @param rNodeGroup node group identifier
//! @param rCoordinates matrix (rows/nodes columns/rCoordinates)
void NuTo::StructureBase::NodeGroupGetCoordinates(int rGroupIdent, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoordinates)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::NodeGroupGetCoordinates] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
    	throw MechanicsException("[NuTo::StructureBase::NodeGroupGetCoordinates] Group is not a node group.");
    Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();
    assert(nodeGroup!=0);

    //all nodes have to have the same dimension
    if(nodeGroup->GetNumMembers()<1)
    	throw MechanicsException("[NuTo::StructureBase::NodeGroupGetCoordinates] Group has no members.");

    int numCoords= nodeGroup->begin()->second->GetNumCoordinates();
    //resize the matrix
    rCoordinates.Resize(nodeGroup->GetNumMembers(),numCoords);
	double coord[3];
	int theNode(0);
    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++, theNode++)
    {
		try
		{
			switch (numCoords)
			{
			case 1:
				itNode->second->GetCoordinates1D(coord);
				rCoordinates(theNode,0)=coord[0];
			break;
			case 2:
				itNode->second->GetCoordinates2D(coord);
				rCoordinates(theNode,0)=coord[0];
				rCoordinates(theNode,1)=coord[1];
			break;
			case 3:
				itNode->second->GetCoordinates3D(coord);
				rCoordinates(theNode,0)=coord[0];
				rCoordinates(theNode,1)=coord[1];
				rCoordinates(theNode,2)=coord[2];
			break;
			default:
				throw MechanicsException("[NuTo::StructureBase::NodeGroupGetCoordinates] The number of coordinates components is either 1, 2 or 3.");
			}
		}
		catch(NuTo::MechanicsException & b)
		{
			b.AddMessage("[NuTo::StructureBase::NodeGroupGetCoordinates] Error getting coordinates.");
			throw b;
		}
		catch(...)
		{
			throw MechanicsException("[NuTo::StructureBase::NodeGroupGetCoordinates] Error getting coordinates of node (unspecified exception).");
		}
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::NodeGroupGetCoordinates] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief gets the global nonlocal eq plastic strain variables of a node
//! @param rNode node identifier
//! @return global (nodal) nonlocal eq plastic strain
void NuTo::StructureBase::NodeGetNonlocalEqPlasticStrain(int rNode, NuTo::FullVector<double,Eigen::Dynamic>& rNonlocalEqPlasticStrain)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	const NodeBase* nodePtr = NodeGetNodePtr(rNode);

	try
	{
		if (nodePtr->GetNumNonlocalEqPlasticStrain()!=2)
		{
			throw MechanicsException("[NuTo::StructureBase::NodeGetNonlocalEqPlasticStrain] Node does not have nonlocal equivalent plastic strains.");
		}
		rNonlocalEqPlasticStrain.resize(2);
		nodePtr->GetNonlocalEqPlasticStrain(rNonlocalEqPlasticStrain.data());
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeGetNonlocalEqPlasticStrain] Error getting global damage.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeGetNonlocalEqPlasticStrain] Error getting NodeGetNonlocalEqPlasticStrain of node (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureBase::NodeGetNonlocalEqPlasticStrain] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief gets the global nonlocal total strain variables of a node
//! @param rNode node identifier
//! @return global (nodal) nonlocal total strain
void NuTo::StructureBase::NodeGetNonlocalTotalStrain(int rNode, NuTo::FullVector<double,Eigen::Dynamic>& rNonlocalTotalStrain)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	const NodeBase* nodePtr = NodeGetNodePtr(rNode);

	try
	{
		switch (nodePtr->GetNumNonlocalTotalStrain())
		{
		case 0:
			throw MechanicsException("[NuTo::StructureBase::NodeGetNonlocalTotalStrain] Node does not have nonlocal total strains.");
        break;
		case 1:
			rNonlocalTotalStrain.resize(1);
			nodePtr->GetNonlocalTotalStrain1D(rNonlocalTotalStrain.data());
		break;
		case 3:
			rNonlocalTotalStrain.resize(3);
			nodePtr->GetNonlocalTotalStrain2D(rNonlocalTotalStrain.data());
		break;
		case 6:
			rNonlocalTotalStrain.resize(6);
			nodePtr->GetNonlocalTotalStrain3D(rNonlocalTotalStrain.data());
		break;
		default:
			throw MechanicsException("[NuTo::StructureBase::NodeGetNonlocalTotalStrain] Number of nonlocal total strain components is either 1, 3 or 6 .");
		break;
		}
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeGetNonlocalTotalStrain] Error getting nonlocal total strain.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeGetNonlocalTotalStrain] Error getting nonlocal total strain of node (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureBase::NodeGetNonlocalTotalStrain] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief calculate the internal force vector for a node
//! @param rId ... node id
//! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
void NuTo::StructureBase::NodeInternalForce(int rId, NuTo::FullVector<double,Eigen::Dynamic>& rNodeForce)
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
void NuTo::StructureBase::NodeGroupInternalForce(int rGroupIdent, NuTo::FullVector<double,Eigen::Dynamic>& rNodeForce)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::NodeGroupInternalForce] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
    	throw MechanicsException("[NuTo::StructureBase::NodeGroupInternalForce] Group is not a node group.");
    const Group<NodeBase> *nodeGroup = dynamic_cast<const Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

	NuTo::FullVector<double,Eigen::Dynamic>nodeForceLocal;

	if (nodeGroup->GetNumMembers()==0)
		throw MechanicsException("[NuTo::StructureBase::NodeGroupInternalForce] Node group is empty.");
	rNodeForce.Resize(nodeGroup->begin()->second->GetNumDisplacements());

    for (Group<NodeBase>::const_iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
		try
		{
			NodeInternalForce(itNode->second, nodeForceLocal);
			if (nodeForceLocal.GetNumRows()!=rNodeForce.GetNumRows())
				throw MechanicsException("[NuTo::StructureBase::NodeGroupInternalForce] The number of displacement components is not equal for all members of the group.");
			rNodeForce+=nodeForceLocal;
		}
		catch(NuTo::MechanicsException & b)
		{
			b.AddMessage("[NuTo::StructureBase::NodeGroupInternalForce] Error getting gradient of internal potential.");
			throw b;
		}
		catch(...)
		{
			throw MechanicsException("[NuTo::StructureBase::NodeGroupInternalForce] Error getting gradient of internal potential (unspecified exception).");
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
void NuTo::StructureBase::NodeInternalForce(const NodeBase* rNodePtr, NuTo::FullVector<double,Eigen::Dynamic>& rNodeForce)
{
	try
	{
		boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;

		boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>( elementOutput )( Element::INTERNAL_GRADIENT );
		boost::assign::ptr_map_insert<ElementOutputVectorInt>( elementOutput )( Element::GLOBAL_ROW_DOF );

		rNodeForce.Resize(rNodePtr->GetNumDisplacements());

		//go through all elements and check, if the node belongs to the element
		std::vector<ElementBase*> elements;
		GetElementsTotal(elements);
		for (unsigned int countElement=0; countElement<elements.size(); countElement++)
		{
			ElementBase* elementPtr=elements[countElement];
			for (int countNode=0; countNode<elementPtr->GetNumNodes(); countNode++)
			{
				if (elementPtr->GetNode(countNode)==rNodePtr)
				{
					elementPtr->Evaluate(elementOutput);

					NuTo::FullVector<double,Eigen::Dynamic>&  elementVector(elementOutput.find(Element::INTERNAL_GRADIENT)->second->GetFullVectorDouble());
	    			std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

	    			assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

					for (int countDof=0; countDof< rNodePtr->GetNumDisplacements(); countDof++)
					{
                        int theDof = rNodePtr->GetDofDisplacement(countDof);
                        for (unsigned int countGlobalDofs=0; countGlobalDofs<elementVectorGlobalDofs.size(); countGlobalDofs++)
                        {
                        	if (elementVectorGlobalDofs[countGlobalDofs] == theDof)
                        	{
                        		rNodeForce(countDof)+=elementVector(countGlobalDofs);
                        	}
                        }
					}
				}
			}
		}
	}
    catch(NuTo::MechanicsException & b)
	{
        b.AddMessage("[NuTo::StructureBase::NodeInternalForce] Error getting gradient of internal potential.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureBase::NodeInternalForce] Error getting gradient of internal potential (unspecified exception).");
	}

}

//! @brief ... store all element ids connected to this node in a vector
//! @param rNode (Input) 			... node id
//! @param rElementNumbers (Output) ... vector of element ids
void NuTo::StructureBase::NodeGetElements(const int rNodeId, NuTo::FullVector<int,Eigen::Dynamic>& rElementNumbers)
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

#ifdef ENABLE_VISUALIZE
//! @brief ... adds all the nodes in the vector to the data structure that is finally visualized
void NuTo::StructureBase::NodeTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{
    std::vector<const NodeBase*> nodeVec;
    this->GetNodesTotal(nodeVec);
    NodeVectorAddToVisualize(rVisualize,rWhat,nodeVec);
}

//! @brief ... adds all the nodes in the vector to the data structure that is finally visualized
void NuTo::StructureBase::NodeVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, const std::vector<const NodeBase*>& rNodes) const
{
    for (unsigned int nodeCount = 0; nodeCount < rNodes.size(); nodeCount++)
    {
        rNodes[nodeCount]->Visualize(rVisualize, rWhat);
    }
}
#endif //ENABLE_VISUALIZE
