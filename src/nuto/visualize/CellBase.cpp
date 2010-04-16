// $Id$

#include "nuto/visualize/CellBase.h"
#include "nuto/visualize/VisualizeException.h"
#include "nuto/visualize/VisualizeDataScalar.h"
#include "nuto/visualize/VisualizeDataVector.h"
#include "nuto/visualize/VisualizeDataTensor.h"
#include "nuto/visualize/VisualizeDataField.h"

// constructor
NuTo::CellBase::CellBase(const std::vector<VisualizeDataType>& rDataTypes)
{
    for (unsigned int DataCount = 0; DataCount < rDataTypes.size(); DataCount++)
    {
        switch (rDataTypes[DataCount].GetDataType())
        {
        case NuTo::VisualizeDataType::SCALAR:
            this->AddDataScalar(DataCount);
            break;
        case NuTo::VisualizeDataType::VECTOR:
            this->AddDataVector(DataCount);
            break;
        case NuTo::VisualizeDataType::TENSOR:
            this->AddDataTensor(DataCount);
            break;
        case NuTo::VisualizeDataType::FIELD:
            this->AddDataField(DataCount, rDataTypes[DataCount].GetNumData());
            break;
        default:
            throw NuTo::VisualizeException("[NuTo::CellBase::CellBase] Unsupported data type.");
        }
    }
}

// get point data
const NuTo::VisualizeDataBase* NuTo::CellBase::GetData(unsigned int rDataIndex) const
{
    if (rDataIndex >= this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::CellBase::GetData] invalid data index.");
    }
    return &(this->mData[rDataIndex]);
}

// add scalar data
void NuTo::CellBase::AddDataScalar(unsigned int rDataIndex)
{
    if (rDataIndex != this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::CellBase::AddDataScalar] invalid data index.");
    }
    this->mData.push_back(new NuTo::VisualizeDataScalar());
}

// add vector data
void NuTo::CellBase::AddDataVector(unsigned int rDataIndex)
{
    if (rDataIndex != this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::CellBase::AddDataVector] invalid data index.");
    }
    this->mData.push_back(new NuTo::VisualizeDataVector());
}

// add tensor data
void NuTo::CellBase::AddDataTensor(unsigned int rDataIndex)
{
    if (rDataIndex != this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::CellBase::AddDataTensor] invalid data index.");
    }
    this->mData.push_back(new NuTo::VisualizeDataTensor());
}

// set tensor data
void NuTo::CellBase::SetDataTensor(unsigned int rDataIndex, double rData[9])
{
    if (rDataIndex >= this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::CellBase::SetDataTensor] invalid data index.");
    }
    if (this->mData[rDataIndex].GetDataType() != NuTo::VisualizeDataType::TENSOR)
    {
        throw NuTo::VisualizeException("[NuTo::CellBase::SetDataTensor] invalid data type.");
    }
    this->mData[rDataIndex].SetData(rData);
}

// add field data
void NuTo::CellBase::AddDataField(unsigned int rDataIndex, unsigned int rNumData)
{
    if (rDataIndex != this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::CellBase::AddDataField] invalid data index.");
    }
    this->mData.push_back(new NuTo::VisualizeDataField(rNumData));
}
