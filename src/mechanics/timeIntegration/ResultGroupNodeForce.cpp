/*
 * ResultDispNode.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */

#include "mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

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

    return groupNode.begin()->second->GetNum(Node::eDof::DISPLACEMENTS);
}

NuTo::eTimeIntegrationResultType NuTo::ResultGroupNodeForce::GetResultType() const
{
    return NuTo::eTimeIntegrationResultType::GROUP_NODE_FORCE;
}


void NuTo::ResultGroupNodeForce::CalculateValues(const StructureBase& rStructure, const Eigen::VectorXd& rResidual_j, const Eigen::VectorXd& rResidual_k, Eigen::MatrixXd& rResult) const
{
    const Group<NodeBase>& nodeGroup = *rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode();

    rResult.resize(1, rStructure.GetDimension());

    for (auto& itNode : nodeGroup)
    {
        assert(itNode.second->GetNum(Node::eDof::DISPLACEMENTS) == rStructure.GetDimension());

        for (int iDim = 0; iDim < rStructure.GetDimension(); iDim++)
        {
            int theDof = itNode.second->GetDof(Node::eDof::DISPLACEMENTS, iDim);

            if (theDof < rStructure.GetNumActiveDofs(Node::eDof::DISPLACEMENTS))
            {
                rResult(0, iDim) += rResidual_j(theDof);
            }
            else
            {
                rResult(0, iDim) += rResidual_k(theDof - rStructure.GetNumActiveDofs(Node::eDof::DISPLACEMENTS));
            }
        }
    }
}