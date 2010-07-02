// $Id: ConstitutiveTangentLocal6x6.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#include <cstring>
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal6x6.h"

// constructor
NuTo::ConstitutiveTangentLocal6x6::ConstitutiveTangentLocal6x6() : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
	for (int count=0; count<36; count++)
		mTangent[count]=0.;
}

// destructor
NuTo::ConstitutiveTangentLocal6x6::~ConstitutiveTangentLocal6x6()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentLocal6x6::GetNumberOfRows() const
{
    return 6;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentLocal6x6::GetNumberOfColumns() const
{
    return 6;
}

// get tangent
const double* NuTo::ConstitutiveTangentLocal6x6::GetData() const
{
    return mTangent;
}

