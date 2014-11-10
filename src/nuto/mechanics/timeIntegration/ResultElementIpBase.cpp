/*
 * ResultElementIpBase.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */
#include "nuto/mechanics/timeIntegration/ResultElementIpBase.h"

NuTo::ResultElementIpBase::ResultElementIpBase(const std::string& rIdent, int rElementId) : ResultBase(rIdent)
{
	mElementId = rElementId;
}

void NuTo::ResultElementIpBase::Info() const
{
}

void NuTo::ResultElementIpBase::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot)
{
	assert(rTimeStepPlot>=0);
	FullMatrix<double,1,Eigen::Dynamic> ipValues(1,this->GetNumData(rStructure));
	this->CalculateValues(rStructure,ipValues);
	if (rTimeStepPlot>=mData.GetNumRows())
	{
		this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
	}
	if (ipValues.GetNumColumns()!=mData.GetNumColumns())
		throw MechanicsException("[NuTo::ResultElementIpBase::CalculateAndAddValues] the allocated number of columns is wrong.");
	mData.SetRow(rTimeStepPlot,ipValues);
}
