// $Id$
#include <iostream>
#include "nuto/visualize/VisualizeDataType.h"
#include "nuto/visualize/VisualizeException.h"

NuTo::VisualizeDataType::VisualizeDataType(const std::string& rIdent, eVisualizeDataType rDataType) : mIdent(rIdent), mDataType(rDataType)
{
    switch (rDataType)
    {
    case NuTo::eVisualizeDataType::SCALAR:
        this->mNumData = 1;
        break;
    case NuTo::eVisualizeDataType::VECTOR:
        this->mNumData = 3;
        break;
    case NuTo::eVisualizeDataType::TENSOR:
        this->mNumData = 9;
        break;
    default:
        this->mNumData = 0;
    }
}

void NuTo::VisualizeDataType::SetNumData(unsigned int& rNumData)
{
    if (this->mDataType != NuTo::eVisualizeDataType::FIELD)
    {
        throw NuTo::VisualizeException("[NuTo::DataType::SetNumEntries] Number of entries can only be changed for FIELD data.");
    }
    this->mNumData = rNumData;
}
