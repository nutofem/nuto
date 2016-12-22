// $Id$

#include "visualize/VisualizeDataTensor.h"
#include "visualize/VisualizeDataType.h"

// constructor
NuTo::VisualizeDataTensor::VisualizeDataTensor()
{
    for (unsigned int count = 0; count < 9; count++)
    {
        this->mData[count] = 0.0;
    }
}

// get data type
NuTo::eVisualizeDataType NuTo::VisualizeDataTensor::GetDataType() const
{
    return NuTo::eVisualizeDataType::TENSOR;
}

// get number of data
unsigned int NuTo::VisualizeDataTensor::GetNumData() const
{
    return 9;
}

// get data
const double* NuTo::VisualizeDataTensor::GetData() const
{
    return this->mData;
}

// set data
void NuTo::VisualizeDataTensor::SetData(const double* rData)
{
    for (unsigned int count = 0; count < 9; count++)
    {
        this->mData[count] = rData[count];
    }
}

// output stream
std::ostream& NuTo::VisualizeDataTensor::Output(std::ostream& os) const
{
    os << this->mData[0];
    for (unsigned int count = 1; count < 9; count++)
    {
        os << " " << this->mData[count];
    }
    return os;
}

