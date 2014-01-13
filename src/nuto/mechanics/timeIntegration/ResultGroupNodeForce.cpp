/*
 * ResultDispNode.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */

#include "nuto/mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"

NuTo::ResultGroupNodeForce::ResultGroupNodeForce(const std::string& rIdent, int rGroupNodeId) : ResultGroupNodeDof(rIdent, rGroupNodeId)
{
}

//! @brief number of dofs (e.g. number of displacement components of a node
int NuTo::ResultGroupNodeForce::GetNumData(const StructureBase& rStructure)const
{
	const Group<NodeBase>* groupNode(rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode());

    //all nodes have to have the same dimension (number of displacement components)
    if(groupNode->GetNumMembers()<1)
    	throw MechanicsException("[NuTo::ResultGroupNodeForce::GetNumData] Group has no members.");

	return groupNode->begin()->second->GetNumDisplacements();
}


void NuTo::ResultGroupNodeForce::CalculateValues(const StructureBase& rStructure,
		const FullVector<double,Eigen::Dynamic>& rResidual_j,
		const FullVector<double,Eigen::Dynamic>& rResidual_k,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)const
{
	const Group<NodeBase>* nodeGroup = rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode();

	rResult.Resize(1,rStructure.GetDimension());

	for (Group<NodeBase>::const_iterator itNode = nodeGroup->begin();
			itNode != nodeGroup->end(); itNode++)
	{
		assert(itNode->second->GetNumDisplacements()==rStructure.GetDimension());

		for (int countDimension = 0; countDimension < rStructure.GetDimension(); countDimension++)
		{
			int theDof = itNode->second->GetDofDisplacement(
					countDimension);
			if (theDof < rStructure.GetNumActiveDofs())
			{
				rResult(0, countDimension) += rResidual_j(
						theDof);
			}
			else
			{
				rResult(0, countDimension) += rResidual_k(
						theDof - rStructure.GetNumActiveDofs());
			}
		}
	}
}
