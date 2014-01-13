
#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/filesystem.hpp>
#include "nuto/mechanics/timeIntegration/ResultBase.h"
#include "nuto/mechanics/structures/StructureBase.h"
//! @brief constructor
NuTo::ResultBase::ResultBase(const std::string& rIdent)
{
	mIdent = rIdent;
}

//! @brief deconstructor
NuTo::ResultBase::~ResultBase()
{
}

void NuTo::ResultBase::SetIdent(const std::string& rIdent)
{
	mIdent = rIdent;
}

std::string NuTo::ResultBase::GetIdent()const
{
	return mIdent;
}

void NuTo::ResultBase::WriteToFile(const std::string& rResultDir, int rTimeStepPlot)const
{
	boost::filesystem::path resultFileName(rResultDir);
	resultFileName /= mIdent+".dat";
	mData.GetBlock(0,0,rTimeStepPlot+1,mData.GetNumColumns()).WriteToFile(resultFileName.string(), "  ");
}

void NuTo::ResultBase::Resize(const StructureBase& rStructure, int rNumTimeSteps, bool rInitValues)
{
	if (rInitValues==true)
	{
		mData.Resize(rNumTimeSteps,this->GetNumData(rStructure));
	}
	else
	{
		mData.ConservativeResize(rNumTimeSteps,this->GetNumData(rStructure));
	}
}
