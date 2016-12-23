
#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/filesystem.hpp>
#include "nuto/math/FullMatrix.h"
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

    int numRowsPerStep = mData.GetNumRows()/(rTimeStepPlot+1); // outputted rows per time step
    mData.GetBlock(numRowsPerStep*rTimeStepPlot, 0, numRowsPerStep, mData.GetNumColumns()).AppendToFile(resultFileName.string(), "  ");
}

void NuTo::ResultBase::Resize(const StructureBase& rStructure, int rNumTimeSteps, bool rInitValues)
{
    int rows(0), cols(0);

    this->GetNumData(rStructure, rows, cols);

	if (rInitValues==true)
	{
        mData.Resize(rNumTimeSteps*rows, cols);
	}
	else
	{
        mData.ConservativeResize(rNumTimeSteps*rows, cols);
	}
}
