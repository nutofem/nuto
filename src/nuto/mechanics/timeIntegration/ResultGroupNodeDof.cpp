/*
 * ResultGroupNodeDof.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/ResultGroupNodeDof.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"

NuTo::ResultGroupNodeDof::ResultGroupNodeDof(const std::string& rIdent, int rGroupNodeId) : ResultBase(rIdent)
{
	mGroupNodeId = rGroupNodeId;
}

void NuTo::ResultGroupNodeDof::Info() const
{
}

void NuTo::ResultGroupNodeDof::CalculateAndAddValues(const NuTo::StructureBase& rStructure, int rTimeStepPlot,
		const FullVector<double,Eigen::Dynamic>& rResidual_j,
		const FullVector<double,Eigen::Dynamic>& rResidual_k)
{
	assert(rTimeStepPlot>=0);
	if (rTimeStepPlot>=mData.GetNumRows())
	{
		this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
	}
	FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> values(1,this->GetNumData(rStructure));
	CalculateValues(rStructure,rResidual_j,rResidual_k, values);

	if (values.GetNumColumns()!=mData.GetNumColumns())
		throw MechanicsException("[NuTo::ResultGroupNodeDof::CalculateAndAddValues] the allocated number of columns is wrong.");

	mData.SetRow(rTimeStepPlot,values);

}

const NuTo::Group<NuTo::NodeBase>* NuTo::ResultGroupNodeDof::GetGroupNodePtr(const StructureBase& rStructure)const
{
	return rStructure.GroupGetGroupPtr(mGroupNodeId)->AsGroupNode();
}
