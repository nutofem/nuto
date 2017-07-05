
#include <boost/filesystem.hpp>
#include "math/EigenCompanion.h"

#include "mechanics/timeIntegration/ResultBase.h"
#include "mechanics/structures/StructureBase.h"

NuTo::ResultBase::ResultBase(const std::string& rIdent)
{
    mIdent = rIdent;
}

NuTo::ResultBase::~ResultBase()
{
}

void NuTo::ResultBase::SetIdent(const std::string& rIdent)
{
    mIdent = rIdent;
}

std::string NuTo::ResultBase::GetIdent() const
{
    return mIdent;
}

void NuTo::ResultBase::WriteToFile(const std::string& rResultDir, int rTimeStepPlot) const
{
    boost::filesystem::path resultFileName(rResultDir);
    resultFileName /= mIdent + ".dat";
    NuTo::EigenCompanion::WriteToFile(mData.block(0, 0, rTimeStepPlot + 1, mData.cols()), resultFileName.string(),
                                      "  ");
}

void NuTo::ResultBase::Resize(const StructureBase& rStructure, int rNumTimeSteps, bool rInitValues)
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
