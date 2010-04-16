// $Id$

#include "nuto/visualize/VisualizeDataType.h"
#include "nuto/visualize/VisualizeException.h"

NuTo::VisualizeDataType::VisualizeDataType(const std::string& rIdent, eDataType rDataType) : mIdent(rIdent), mDataType(rDataType)
{
    switch (rDataType)
    {
    case NuTo::VisualizeDataType::SCALAR:
        this->mNumData = 1;
        break;
    case NuTo::VisualizeDataType::VECTOR:
        this->mNumData = 3;
        break;
    case NuTo::VisualizeDataType::TENSOR:
        this->mNumData = 9;
        break;
    default:
        this->mNumData = 0;
    }
}

void NuTo::VisualizeDataType::SetNumData(unsigned int& rNumData)
{
    if (this->mDataType != NuTo::VisualizeDataType::FIELD)
    {
        throw NuTo::VisualizeException("[NuTo::DataType::SetNumEntries] Number of entries can only be changed for FIELD data.");
    }
    this->mNumData = rNumData;
}
