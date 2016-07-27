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
	const Group<NodeBase>& groupNode = *rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode();

    //all nodes have to have the same dimension (number of displacement components)
    if(groupNode.GetNumMembers() < 1)
    	throw MechanicsException("[NuTo::ResultGroupNodeForce::GetNumData] Group has no members.");

	return groupNode.begin()->second->GetNum(Node::DISPLACEMENTS);
}


void NuTo::ResultGroupNodeForce::CalculateValues(const StructureBase& rStructure, const FullVector<double, Eigen::Dynamic>& rResidual_j, const FullVector<double, Eigen::Dynamic>& rResidual_k, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const
{
    const Group<NodeBase>& nodeGroup = *rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode();

    rResult.Resize(1, rStructure.GetDimension());

    for (auto& itNode : nodeGroup)
    {
        assert(itNode.second->GetNum(Node::DISPLACEMENTS) == rStructure.GetDimension());

        for (int iDim = 0; iDim < rStructure.GetDimension(); iDim++)
        {
            int theDof = itNode.second->GetDof(Node::DISPLACEMENTS, iDim);

            if (theDof < rStructure.GetNumActiveDofs(Node::DISPLACEMENTS))
            {
                rResult(0, iDim) += rResidual_j(theDof);
            }
            else
            {
                rResult(0, iDim) += rResidual_k(theDof - rStructure.GetNumActiveDofs(Node::DISPLACEMENTS));
            }
        }
    }
}
