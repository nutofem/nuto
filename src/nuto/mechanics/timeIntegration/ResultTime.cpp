/*
 * ResultTime.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */
#include "nuto/mechanics/timeIntegration/ResultTime.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"

NuTo::ResultTime::ResultTime(const std::string& rIdent) : ResultBase(rIdent)
{
}

void NuTo::ResultTime::Info() const
{
}

void NuTo::ResultTime::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot, double rTime)
{
	assert(rTimeStepPlot>=0);
	if (rTimeStepPlot>=mData.GetNumRows())
	{
		this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
	}
	mData(rTimeStepPlot,0) = rTime;
}

NuTo::eTimeIntegrationResultType NuTo::ResultTime::GetResultType() const
{
    return NuTo::eTimeIntegrationResultType::TIME;
}
