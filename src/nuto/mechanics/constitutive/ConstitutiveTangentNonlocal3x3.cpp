// $Id: ConstitutiveTangentNonlocal3x3.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#include <cstring>
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal3x3.h"

// constructor
NuTo::ConstitutiveTangentNonlocal3x3::ConstitutiveTangentNonlocal3x3(int rNumNonlocalElements) : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
	mNonlocalMatrices.resize(rNumNonlocalElements);
}

// destructor
NuTo::ConstitutiveTangentNonlocal3x3::~ConstitutiveTangentNonlocal3x3()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentNonlocal3x3::GetNumberOfRows() const
{
    return mNonlocalMatrices.size() * 3;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentNonlocal3x3::GetNumberOfColumns() const
{
    return 3;
}

// get tangent
const double* NuTo::ConstitutiveTangentNonlocal3x3::GetData() const
{
    throw MechanicsException("[NuTo::ConstitutiveTangentNonlocal3x3::GetData] Method not implemented for nonlocal matrix.");
}

//! @brief ... get a local submatrix
//! @param ... rSubMatrix number of the submatrix
//! @return ... pointer to the tangent submatrix 3x3 matrix (column major storage)
NuTo::ConstitutiveTangentLocal3x3* NuTo::ConstitutiveTangentNonlocal3x3::GetSubMatrix(int rSubMatrix)
{
	assert(rSubMatrix<(int)mNonlocalMatrices.size());
	return &mNonlocalMatrices[rSubMatrix];
}

//! @brief ... get a local submatrix
//! @param ... rSubMatrix number of the submatrix
//! @return ... pointer to the tangent submatrix 3x3 matrix (column major storage)
const NuTo::ConstitutiveTangentLocal3x3* NuTo::ConstitutiveTangentNonlocal3x3::GetSubMatrix(int rSubMatrix)const
{
	assert(rSubMatrix<(int)mNonlocalMatrices.size());
	return &mNonlocalMatrices[rSubMatrix];
}

//! @brief ... get a the number of local submatrices
//! @return ... number of submatrices
int NuTo::ConstitutiveTangentNonlocal3x3::GetNumSubMatrices()const
{
	return (int)mNonlocalMatrices.size();
}
