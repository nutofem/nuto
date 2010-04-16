// $Id$

#include "nuto/visualize/VisualizeDataScalar.h"

// constructor
NuTo::VisualizeDataScalar::VisualizeDataScalar()
{
    this->mData[0] = 0;
}

// get data type
NuTo::VisualizeDataType::eDataType NuTo::VisualizeDataScalar::GetDataType() const
{
    return NuTo::VisualizeDataType::SCALAR;
}

// get number of data
unsigned int NuTo::VisualizeDataScalar::GetNumData() const
{
    return 1;
}

// get data
const double* NuTo::VisualizeDataScalar::GetData() const
{
    return this->mData;
}

// set data
void NuTo::VisualizeDataScalar::SetData(const double* rData)
{
    this->mData[0] = rData[0];
}

// output stream
std::ostream& NuTo::VisualizeDataScalar::Output(std::ostream& os) const
{
    os << this->mData[0];
    return os;
}

