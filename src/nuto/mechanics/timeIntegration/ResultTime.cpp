/*
 * ResultTime.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */
#include <boost/filesystem.hpp>
#include "nuto/mechanics/timeIntegration/ResultTime.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::ResultTime::ResultTime(const std::string& rIdent) : ResultBase(rIdent)
{
}

void NuTo::ResultTime::Info() const
{
}

void NuTo::ResultTime::Resize(int rNumTimeSteps, bool rInitValues)
{
	if (rInitValues==true)
	{
		mData.Resize(rNumTimeSteps);
	}
	else
	{
		mData.ConservativeResize(rNumTimeSteps);
	}
}

void NuTo::ResultTime::CalculateAndAddValues(int rTimeStepPlot, double rTime)
{
	assert(rTimeStepPlot>=0);
	if (rTimeStepPlot>=mData.GetNumRows())
	{
		this->Resize(2*(rTimeStepPlot+1),false);
	}
	mData(rTimeStepPlot) = rTime;
}

void NuTo::ResultTime::WriteToFile(const std::string& rResultDir, int rTimeStepPlot)const
{
	boost::filesystem::path resultFileName(rResultDir);
	resultFileName /= mIdent+".dat";
	mData.GetBlock(0,0,rTimeStepPlot+1,1).WriteToFile(resultFileName.string(), "  ");
}
