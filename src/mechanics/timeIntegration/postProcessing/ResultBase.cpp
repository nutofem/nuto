#include "mechanics/timeIntegration/postProcessing/ResultBase.h"

#include <boost/filesystem.hpp>
#include "math/EigenCompanion.h"

#include "mechanics/structures/StructureBase.h"

using namespace NuTo;

ResultBase::ResultBase(const std::string& rIdent)
{
    mIdent = rIdent;
}


ResultBase::~ResultBase()
{
}


void ResultBase::SetIdent(const std::string& rIdent)
{
    mIdent = rIdent;
}


std::string ResultBase::GetIdent() const
{
    return mIdent;
}


void ResultBase::WriteToFile(const std::string& rResultDir, int rTimeStepPlot) const
{
    boost::filesystem::path resultFileName(rResultDir);
    resultFileName /= mIdent + ".dat";
    EigenCompanion::WriteToFile(mData.block(0, 0, rTimeStepPlot + 1, mData.cols()), resultFileName.string(), "  ");
}


void ResultBase::Resize(const StructureBase& rStructure, int rNumTimeSteps, bool rInitValues)
{
    if (rInitValues)
    {
        mData.setZero(rNumTimeSteps, this->GetNumData(rStructure));
    }
    else
    {
        mData.conservativeResize(rNumTimeSteps, this->GetNumData(rStructure));
    }
}

ResultBase* new_clone(const ResultBase& result)
{
    return result.Clone().release();
}
