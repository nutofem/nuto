// $Id$

#include "visualize/VisualizeDataField.h"
#include "visualize/VisualizeDataType.h"

// constructor
NuTo::VisualizeDataField::VisualizeDataField(unsigned int rNumData) : mNumData(rNumData)
{
    if (rNumData > 0)
    {
        this->mData = new double[rNumData];
        for (unsigned int count = 0; count < rNumData; count++)
        {
            this->mData[count] = 0.0;
        }
    }
    else
    {
        this->mData = 0;
    }
}

// destructor
NuTo::VisualizeDataField::~VisualizeDataField()
{
    if (this->mData != 0)
    {
        delete[] this->mData;
    }
}

// get data type
NuTo::eVisualizeDataType NuTo::VisualizeDataField::GetDataType() const
{
    return NuTo::eVisualizeDataType::FIELD;
}

// get number of data
unsigned int NuTo::VisualizeDataField::GetNumData() const
{
    return this->mNumData;
}

// get data
const double* NuTo::VisualizeDataField::GetData() const
{
    return this->mData;
}

// set data
void NuTo::VisualizeDataField::SetData(const double* rData)
{
    for (unsigned int count = 0; count < this->mNumData; count++)
    {
        this->mData[count] = rData[count];
    }
}

// output stream
std::ostream& NuTo::VisualizeDataField::Output(std::ostream& os) const
{
    if (this->mNumData > 0)
    {
        os << this->mData[0];
        for (unsigned int count = 1; count < this->mNumData; count++)
        {
            os << " " << this->mData[count];
        }
    }
    return os;
}

