/*
 * ResultGroupNodeDof.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */

#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/ResultGroupNodeDof.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"

NuTo::ResultGroupNodeDof::ResultGroupNodeDof(const std::string& rIdent, int rGroupNodeId) : ResultBase(rIdent)
{
	mGroupNodeId = rGroupNodeId;
}

void NuTo::ResultGroupNodeDof::Info() const
{
}

void NuTo::ResultGroupNodeDof::CalculateAndAddValues(const NuTo::StructureBase& rStructure, int rTimeStepPlot,
		const Eigen::VectorXd& rResidual_j,
		const Eigen::VectorXd& rResidual_k)
{
	assert(rTimeStepPlot>=0);
	if (rTimeStepPlot>=mData.rows())
	{
		this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
	}
	Eigen::VectorXd values = CalculateValues(rStructure,rResidual_j,rResidual_k);

	if (values.rows()!=mData.cols())
		throw MechanicsException(__PRETTY_FUNCTION__, "the allocated number of rows is wrong.");

	mData.row(rTimeStepPlot) = values.transpose();
}

const NuTo::Group<NuTo::NodeBase>* NuTo::ResultGroupNodeDof::GetGroupNodePtr(const StructureBase& rStructure)const
{
	return rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode();
}
