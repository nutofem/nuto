// $Id$

#include "nuto/visualize/Point.h"
#include "nuto/visualize/VisualizeException.h"
#include "nuto/visualize/VisualizeDataScalar.h"
#include "nuto/visualize/VisualizeDataVector.h"
#include "nuto/visualize/VisualizeDataTensor.h"
#include "nuto/visualize/VisualizeDataField.h"

// constructor
NuTo::Point::Point(const double* rCoordinates)
{
    this->mCoordinates[0] = rCoordinates[0];
    this->mCoordinates[1] = rCoordinates[1];
    this->mCoordinates[2] = rCoordinates[2];
}

// add scalar data
void NuTo::Point::AddDataScalar(unsigned int rDataIndex)
{
    if (rDataIndex != this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::Point::AddDataScalar] invalid data index.");
    }
    this->mData.push_back(new NuTo::VisualizeDataScalar());
}

// add vector data
void NuTo::Point::AddDataVector(unsigned int rDataIndex)
{
    if (rDataIndex != this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::Point::AddDataVector] invalid data index.");
    }
    this->mData.push_back(new NuTo::VisualizeDataVector());
}

// set scalar data
void NuTo::Point::SetDataScalar(unsigned int rDataIndex, double rData)
{
    if (rDataIndex >= this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::Point::SetDataVector] invalid data index.");
    }
    if (this->mData[rDataIndex].GetDataType() != NuTo::VisualizeDataType::SCALAR)
    {
        throw NuTo::VisualizeException("[NuTo::Point::SetDataScalar] invalid data type.");
    }
    this->mData[rDataIndex].SetData(&rData);
}

// set vector data
void NuTo::Point::SetDataVector(unsigned int rDataIndex, double rData[3])
{
    if (rDataIndex >= this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::Point::SetDataVector] invalid data index.");
    }
    if (this->mData[rDataIndex].GetDataType() != NuTo::VisualizeDataType::VECTOR)
    {
        throw NuTo::VisualizeException("[NuTo::Point::SetDataVector] invalid data type.");
    }
    this->mData[rDataIndex].SetData(rData);
}


// add tensor data
void NuTo::Point::AddDataTensor(unsigned int rDataIndex)
{
    if (rDataIndex != this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::Point::AddDataTensor] invalid data index.");
    }
    this->mData.push_back(new NuTo::VisualizeDataTensor());
}

// add field data
void NuTo::Point::AddDataField(unsigned int rDataIndex, unsigned int rNumData)
{
    if (rDataIndex != this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::Point::AddDataField] invalid data index.");
    }
    this->mData.push_back(new NuTo::VisualizeDataField(rNumData));
}

// get point data
const NuTo::VisualizeDataBase* NuTo::Point::GetData(unsigned int rDataIndex) const
{
    if (rDataIndex >= this->mData.size())
    {
        throw NuTo::VisualizeException("[NuTo::Point::GetData] invalid data index.");
    }
    return &(this->mData[rDataIndex]);
}
