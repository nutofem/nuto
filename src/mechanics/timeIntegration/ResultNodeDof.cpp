/*
 * ResultNodeDof.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */


#include "mechanics/timeIntegration/ResultNodeDof.h"
#include "mechanics/nodes/NodeBase.h"

NuTo::ResultNodeDof::ResultNodeDof(const std::string& rIdent, int rNodeId)
    : ResultBase(rIdent)
{
    mNodeId = rNodeId;
}

void NuTo::ResultNodeDof::Info() const
{
}

void NuTo::ResultNodeDof::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot)
{
    assert(rTimeStepPlot >= 0);
    Eigen::VectorXd dofValues = this->CalculateValues(rStructure);
	if (rTimeStepPlot>=mData.rows())
	{
		this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
	}
	if (dofValues.rows()!=mData.cols())
		throw Exception(__PRETTY_FUNCTION__, "the allocated number of rows is wrong.");
	mData.row(rTimeStepPlot) = dofValues.transpose();
}
