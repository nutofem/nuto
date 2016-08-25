// $Id$

#include "nuto/visualize/VisualizeDataVector.h"
#include "nuto/visualize/VisualizeDataType.h"

// constructor
NuTo::VisualizeDataVector::VisualizeDataVector()
{
    this->mData[0] = 0.0;
    this->mData[1] = 0.0;
    this->mData[2] = 0.0;
}

// get data type
NuTo::eVisualizeDataType NuTo::VisualizeDataVector::GetDataType() const
{
    return NuTo::eVisualizeDataType::VECTOR;
}

// get number of data
unsigned int NuTo::VisualizeDataVector::GetNumData() const
{
    return 3;
}

// get data
const double* NuTo::VisualizeDataVector::GetData() const
{
    return this->mData;
}

// set data
void NuTo::VisualizeDataVector::SetData(const double* rData)
{
    this->mData[0] = rData[0];
    this->mData[1] = rData[1];
    this->mData[2] = rData[2];
}

// output stream
std::ostream& NuTo::VisualizeDataVector::Output(std::ostream& os) const
{
    os << this->mData[0] << " " << this->mData[1] << " " << this->mData[2];
    return os;
}

