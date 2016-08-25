/*
 * ResultNodeDof.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/timeIntegration/ResultNodeDof.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::ResultNodeDof::ResultNodeDof(const std::string& rIdent, int rNodeId) : ResultBase(rIdent)
{
	mNodeId = rNodeId;
}

void NuTo::ResultNodeDof::Info() const
{
}

void NuTo::ResultNodeDof::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot)
{
	assert(rTimeStepPlot>=0);
	FullMatrix<double,1,Eigen::Dynamic> dofValues(1,this->GetNumData(rStructure));
	this->CalculateValues(rStructure,dofValues);
	if (rTimeStepPlot>=mData.GetNumRows())
	{
		this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
	}
	if (dofValues.GetNumColumns()!=mData.GetNumColumns())
		throw MechanicsException("[NuTo::ResultNodeDof::CalculateAndAddValues] the allocated number of columns is wrong.");
	mData.SetRow(rTimeStepPlot,dofValues);
}
