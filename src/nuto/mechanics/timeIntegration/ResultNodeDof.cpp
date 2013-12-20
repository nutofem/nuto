/*
 * ResultNodeDof.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */
#include <boost/filesystem.hpp>
#include "nuto/mechanics/timeIntegration/ResultNodeDof.h"
#include "nuto/mechanics/nodes/NodeBase.h"

NuTo::ResultNodeDof::ResultNodeDof(const std::string& rIdent, const NodeBase* rNodePtr) : ResultBase(rIdent)
{
	mNodePtr = rNodePtr;
}

void NuTo::ResultNodeDof::Info() const
{
}

void NuTo::ResultNodeDof::Resize(int rNumTimeSteps, bool rInitValues)
{
	if (rInitValues==true)
	{
		mData.Resize(rNumTimeSteps,this->GetNumDofs());
	}
	else
	{
		mData.ConservativeResize(rNumTimeSteps,this->GetNumDofs());
	}
}

void NuTo::ResultNodeDof::CalculateAndAddValues(int rTimeStepPlot)
{
	assert(rTimeStepPlot>=0);
	FullMatrix<double,1,Eigen::Dynamic> dofValues(1,this->GetNumDofs());
	this->CalculateValues(dofValues);
	if (rTimeStepPlot>=mData.GetNumRows())
	{
		this->Resize(2*(rTimeStepPlot+1),false);
	}
	if (dofValues.GetNumColumns()!=mData.GetNumColumns())
		throw MechanicsException("[NuTo::ResultNodeDof::CalculateAndAddValues] the allocated number of columns is wrong.");
	mData.SetRow(rTimeStepPlot,dofValues);
}

void NuTo::ResultNodeDof::WriteToFile(const std::string& rResultDir, int rTimeStepPlot)const
{
	boost::filesystem::path resultFileName(rResultDir);
	resultFileName /= mIdent+".dat";
	mData.GetBlock(0,0,rTimeStepPlot+1,1).WriteToFile(resultFileName.string(), "  ");
}
