// $Id: ConstitutiveTangentLocal3x3.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#include <cstring>
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"

// constructor
NuTo::ConstitutiveTangentLocal3x3::ConstitutiveTangentLocal3x3() : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
    mTangent[0] = 0.;
    mTangent[1] = 0.;
    mTangent[2] = 0.;
    mTangent[3] = 0.;
    mTangent[4] = 0.;
    mTangent[5] = 0.;
    mTangent[6] = 0.;
    mTangent[7] = 0.;
    mTangent[8] = 0.;
}

// destructor
NuTo::ConstitutiveTangentLocal3x3::~ConstitutiveTangentLocal3x3()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentLocal3x3::GetNumberOfRows() const
{
    return 3;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentLocal3x3::GetNumberOfColumns() const
{
    return 3;
}

// get tangent
const double* NuTo::ConstitutiveTangentLocal3x3::GetData() const
{
    return mTangent;
}

