// $Id: ConstitutiveTangentLocal1x1.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#include <cstring>
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"

// constructor
NuTo::ConstitutiveTangentLocal1x1::ConstitutiveTangentLocal1x1() : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
    this->mTangent = 0;
}

// destructor
NuTo::ConstitutiveTangentLocal1x1::~ConstitutiveTangentLocal1x1()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentLocal1x1::GetNumberOfRows() const
{
    return 1;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentLocal1x1::GetNumberOfColumns() const
{
    return 1;
}

// get tangent
const double* NuTo::ConstitutiveTangentLocal1x1::GetData() const
{
    return &mTangent;
}

