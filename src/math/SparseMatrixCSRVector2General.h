// $Id: SparseMatrixCSRVector2General.h 235 2010-04-22 09:25:38Z arnold2 $

#pragma once

#include <algorithm>

#include "math/SparseMatrixCSRVector2General_Def.h"

#include "math/FullMatrix.h"
//#include "math/SparseMatrixCSRGeneral_Def.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRVector2Symmetric.h"
#include "math/MathException.h"

//! @brief ... constructor
//! @param rNumRows_ ... number of rows
//! @param rNumColumns_ ... number of columns
template<class T>
NuTo::SparseMatrixCSRVector2General<T>::SparseMatrixCSRVector2General(int rNumRows, int rNumColumns) :
                             NuTo::SparseMatrixCSRVector2<T>(rNumRows)
{
	mNumColumns = rNumColumns;
}

//! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a predefined tolerance)
//! @param rFullMatrix ... input matrix (full storage)
//! @param rAbsoluteTolerance ... absolute tolerance
//! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance * max(abs(rMatrixEntry))
template<class T>
NuTo::SparseMatrixCSRVector2General<T>::SparseMatrixCSRVector2General(const NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rFullMatrix, double rAbsoluteTolerance, double rRelativeTolerance):
                             NuTo::SparseMatrixCSRVector2<T>(rFullMatrix.GetNumRows())
{
	this->mNumColumns = rFullMatrix.GetNumColumns();

	double tolerance = rAbsoluteTolerance;
	if (rRelativeTolerance > 1e-14)
	{
		T maxValue = this->Max();
		T minValue = this->Min();
		if (std::abs(maxValue)>std::abs(minValue))
			tolerance += rRelativeTolerance * std::abs(maxValue);
		else
			tolerance += rRelativeTolerance * std::abs(minValue);

	}

	for (int row = 0; row < rFullMatrix.GetNumRows(); row++)
	{
		for (int col = 0; col < rFullMatrix.GetNumColumns(); col++)
		{
			if (std::abs(rFullMatrix(row,col)) > tolerance)
			{
				this->AddValue(row,col,rFullMatrix(row,col));
			}
		}
	}
}

//! @brief ... create sparse matrix with vector of vector from standard CSR format
//! @param rCSRMatrix ... input matrix (full storage)
template<class T>
NuTo::SparseMatrixCSRVector2General<T>::SparseMatrixCSRVector2General(const NuTo::SparseMatrixCSRGeneral<T>& rCSRMatrix) :
              NuTo::SparseMatrixCSRVector2<T>::SparseMatrixCSRVector2(rCSRMatrix.GetNumRows())
{
    this->mNumColumns = rCSRMatrix.GetNumColumns();
    this->mColumns.resize(rCSRMatrix.GetNumRows());
    this->mValues.resize(rCSRMatrix.GetNumRows());
    std::vector<int>::const_iterator startIteratorColumns(rCSRMatrix.mColumns.begin());
    typename std::vector<T>::const_iterator startIteratorValues(rCSRMatrix.mValues.begin());
    int numEntries;
    for (unsigned int row = 0; row < this->mColumns.size(); row++)
    {
   	    numEntries = rCSRMatrix.mRowIndex[row+1]-rCSRMatrix.mRowIndex[row];
    	this->mColumns[row].resize(numEntries);
    	this->mValues[row].resize(numEntries);
    	std::copy(startIteratorColumns, startIteratorColumns+numEntries,this->mColumns[row].begin());
    	std::copy(startIteratorValues, startIteratorValues+numEntries,this->mValues[row].begin());
    	startIteratorColumns+=numEntries;
    	startIteratorValues+=numEntries;
     }
}

//! @brief ... create sparse matrix with vector of vector from symmetric sparse matrix
//! @param rCSRVector2Symmetric ... input matrix (full storage)
template<class T>
NuTo::SparseMatrixCSRVector2General<T>::SparseMatrixCSRVector2General(const SparseMatrixCSRVector2Symmetric<T>& rCSRVector2Symmetric) :
            NuTo::SparseMatrixCSRVector2<T>::SparseMatrixCSRVector2(rCSRVector2Symmetric.GetNumRows())
{
    Resize(rCSRVector2Symmetric.GetNumRows(), rCSRVector2Symmetric.GetNumColumns());
    for (int row = 0; row < rCSRVector2Symmetric.GetNumRows(); ++row)
    {
        for (unsigned int pos = 0; pos < rCSRVector2Symmetric.mValues[row].size(); ++pos)
        {
            this->AddValue(row, rCSRVector2Symmetric.mColumns[row][pos], rCSRVector2Symmetric.mValues[row][pos]);
            if (row!=rCSRVector2Symmetric.mColumns[row][pos])
                this->AddValue(rCSRVector2Symmetric.mColumns[row][pos], row, rCSRVector2Symmetric.mValues[row][pos]);
        }
    }
}

//! @brief ... returns whether the matrix is symmetric or unsymmetric
//! @return true if the matrix is symmetric and false if the matrix is unsymmetric
template<class T>
bool NuTo::SparseMatrixCSRVector2General<T>::IsSymmetric() const
{
	return false;
}

//! @brief ... returns the number of columns
//! @return number of columns
template <class T>
int NuTo::SparseMatrixCSRVector2General<T>::GetNumColumns() const
{
	return this->mNumColumns;
}

//! @brief ... add nonzero entry to matrix
//! @param rRow ... row of the nonzero entry (zero based indexing!!!)
//! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
//! @param rValue ... value of the nonzero entry
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::AddValue(int rRow, int rColumn, const T& rValue)
{
    // check bounds via asserts
    assert(rRow < (int)this->mValues.size()     && "row index is out of bounds.");
    assert(rRow >= 0                            && "row index is out of bounds.");

    assert(rColumn < mNumColumns                && "column index is out of bounds.");
    assert(rColumn >= 0                         && "column index is out of bounds.");

    if (this->mOneBasedIndexing)
    {
        rColumn++;
    }
    auto& colVec = this->mColumns[rRow];
    auto& valVec = this->mValues[rRow];

    auto it = lower_bound(colVec.begin(), colVec.end(), rColumn);
    if (it == colVec.end())
    {
        // insert new value at the end of the vectors
        colVec.push_back(rColumn);
        valVec.push_back(rValue);
    }
    else
    {
        auto itValue = valVec.begin() + std::distance(colVec.begin(), it);

        if (*it==rColumn)
        {
            // += to existing value
            *itValue += rValue;
        }
        else
        {
            // insert new value
            colVec.insert(it, rColumn);
            valVec.insert(itValue, rValue);
        }
    }
}

//! @brief ... return the matrix type
template<class T>
NuTo::eSparseMatrixType NuTo::SparseMatrixCSRVector2General<T>::GetSparseMatrixType()const
{
    return NuTo::eSparseMatrixType::CSRVECTOR2GENERAL;
}

//! @brief ... import matrix from slang object stored in  a text file
//! @param rFileName ... file name
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::ImportFromSLangText(const char* rFileName)
{
	throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] to be implemented.");
}

//! @brief ... write nonzero matrix entries into a matrix
//! @param rFullMatrix ... the matrix
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::WriteEntriesToMatrix(Matrix<T>& rMatrix) const
{
	rMatrix.Resize(this->GetNumRows(), this->GetNumColumns());
	if (this->mOneBasedIndexing)
	{
		for (unsigned int row=0; row<this->mColumns.size(); row++)
		{
			for (unsigned int col_count=0; col_count<this->mColumns[row].size(); col_count++)
			{
				rMatrix.AddValue(row, this->mColumns[row][col_count]-1, this->mValues[row][col_count]);
			}
		}
	}
	else
	{
		for (unsigned int row=0; row<this->mColumns.size(); row++)
		{
			for (unsigned int col_count=0; col_count<this->mColumns[row].size(); col_count++)
			{
				rMatrix.AddValue(row, this->mColumns[row][col_count], this->mValues[row][col_count]);
			}
		}
	}
}

//! @brief ... returns the symmetric part of the matrix 0.5*(A+A^T)
//! @return symmetric part
template<class T>
NuTo::SparseMatrixCSRVector2Symmetric<T> NuTo::SparseMatrixCSRVector2General<T>::SymmetricPart() const
{

	if (this->mNumColumns!=(int)this->mValues.size())
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] matrix has to be square.");
	SparseMatrixCSRVector2Symmetric<T> symmetricMatrix(this->mNumColumns,this->mNumColumns);
	if (this->mOneBasedIndexing)
	{
		symmetricMatrix.SetOneBasedIndexing();
		for (unsigned int row=0; row<this->mColumns.size(); row++)
		{
			for (unsigned int col_count=0; col_count<this->mColumns[row].size(); col_count++)
			{
				if (row<=(unsigned int)this->mColumns[row][col_count]-1)
				    symmetricMatrix.AddValue(row, this->mColumns[row][col_count]-1, 0.5*this->mValues[row][col_count]);
				if (row>=(unsigned int)this->mColumns[row][col_count]-1)
				    symmetricMatrix.AddValue(this->mColumns[row][col_count]-1,row, 0.5*this->mValues[row][col_count]);
			}
		}
	}
	else
	{
		symmetricMatrix.SetZeroBasedIndexing();
		for (unsigned int row=0; row<this->mColumns.size(); row++)
		{
			for (unsigned int col_count=0; col_count<this->mColumns[row].size(); col_count++)
			{
				if (row<=(unsigned int)this->mColumns[row][col_count])
				    symmetricMatrix.AddValue(row, this->mColumns[row][col_count], 0.5*this->mValues[row][col_count]);
				if (row>=(unsigned int)this->mColumns[row][col_count])
				    symmetricMatrix.AddValue(this->mColumns[row][col_count],row, 0.5*this->mValues[row][col_count]);
			}
		}
	}
	return symmetricMatrix;
}

//! @brief ... subtract two matrices
//! @param rOther ... general sparse matrix stored in the CSRVector2 format
//! @return reference to this matrix
template<class T>
NuTo::SparseMatrixCSRVector2General<T>& NuTo::SparseMatrixCSRVector2General<T>::operator-=  ( const NuTo::SparseMatrixCSRVector2General<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRVector2General::operator-=] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRVector2General::operator-=] both matrices must have zero based indexing.");
	}
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int row = 0; row < rOther.GetNumRows(); ++row)
    {
        for (unsigned int pos = 0; pos < rOther.mValues[row].size(); ++pos)
        {
            this->AddValue(row, rOther.mColumns[row][pos], -rOther.mValues[row][pos]);
        }
    }
	return *this;
}

//! @brief ... add two matrices
//! @param rOther ... general sparse matrix stored in the CSRVector2 format
//! @return reference to this matrix
template<class T>
NuTo::SparseMatrixCSRVector2General<T>& NuTo::SparseMatrixCSRVector2General<T>::operator+=  ( const NuTo::SparseMatrixCSRVector2General<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");
	}
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int row = 0; row < rOther.GetNumRows(); ++row)
	{
        for (unsigned int pos = 0; pos < rOther.mValues[row].size(); ++pos)
		{
			this->AddValue(row, rOther.mColumns[row][pos], rOther.mValues[row][pos]);
		}
	}
	return *this;
}

//! @brief ... add two matrices
//! @param rOther ... symmetric sparse matrix stored in the CSRVector2 format
//! @return reference to this matrix
template<class T>
NuTo::SparseMatrixCSRVector2General<T>& NuTo::SparseMatrixCSRVector2General<T>::operator+=  ( const NuTo::SparseMatrixCSRVector2Symmetric<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");
	}
    for (int row = 0; row < rOther.GetNumRows(); ++row)
	{
        for (unsigned int pos = 0; pos < rOther.mValues[row].size(); ++pos)
		{
			this->AddValue(row, rOther.mColumns[row][pos], rOther.mValues[row][pos]);
			if (row!=rOther.mColumns[row][pos])
				this->AddValue(rOther.mColumns[row][pos], row, rOther.mValues[row][pos]);
		}
	}
	return *this;
}

//! @brief ... matrix - matrix multiplication
//! @param rOther ... general sparse matrix stored in the CSR format
//! @return general sparse matrix stored in the CSR format
template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::operator* ( const NuTo::SparseMatrixCSRVector2General<T> &rOther ) const
{
	if (this->GetNumColumns() != rOther.GetNumRows())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");
	}
	SparseMatrixCSRVector2General<T> result(this->GetNumRows(), rOther.GetNumColumns());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int thisRow = 0; thisRow < this->GetNumRows(); ++thisRow)
	{
		const std::vector<T>& thisValueVec(this->mValues[thisRow]);
		const std::vector<int>& thisColumnVec(this->mColumns[thisRow]);
        for (unsigned int thisPos = 0; thisPos < thisValueVec.size(); ++thisPos)
		{
			unsigned int thisColumn =thisColumnVec[thisPos];
			T thisValue = thisValueVec[thisPos];

			const std::vector<int>& otherColumnVec(rOther.mColumns[thisColumn]);
			const std::vector<T>& otherValueVec(rOther.mValues[thisColumn]);
            for (unsigned int otherPos = 0; otherPos < otherColumnVec.size(); ++otherPos)
			{
				result.AddValue(thisRow, otherColumnVec[otherPos], thisValue * otherValueVec[otherPos]);
			}
		}
	}
	return result;
}

/*
//! @brief ... multiplies the matrix with an scalar value
//! @param rOther ... scalar value
//! @return ... the multiplied matrix (sparse csr storage)
template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::operator* ( const T &rOther ) const
{
	SparseMatrixCSRVector2General<T> result(*this);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned int thisRow=0; thisRow< this->mValues.size(); ++thisRow)
	{
		std::vector<T>& resultValueVec(result.mValues[thisRow]);
        for (unsigned int thisColumn=0; thisColumn< this->mValues[thisRow].size(); ++thisColumn)
		{
			resultValueVec[thisColumn]*=rOther;
		}
	}
	return result;
}
*/

//! @brief ... multiply sparse matrix with a full matrix
//! @param rFullMatrix ... full matrix which is multiplied with the sparse matrix
//! @return ... full matrix
template<class T>
NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> NuTo::SparseMatrixCSRVector2General<T>::operator* (const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &rMatrix) const
{
	if (this->GetNumColumns() != rMatrix.GetNumRows())
	{
		std::cout << "this->GetNumColumns() " << this->GetNumColumns() << " rMatrix.GetNumRows() " << rMatrix.GetNumRows() << "\n";
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> result(this->GetNumRows(),rMatrix.GetNumColumns());
	if (this->HasOneBasedIndexing())
	{
		// loop over rows
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int row = 0; row < this->GetNumRows(); ++row)
		{
			// initialize result
            for (int matrixCol = 0; matrixCol < result.GetNumColumns(); ++matrixCol)
			{
				result(row,matrixCol) = 0.0;
			}
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			// perform multiplication
            for (unsigned int pos = 0; pos < this->mColumns.size(); ++pos)
			{
				int column = thisColumnVec[pos] - 1;
				T value = thisValueVec[pos];
                for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); ++matrixCol)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	else
	{
		// loop over rows
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int row = 0; row < this->GetNumRows(); ++row)
		{
			// initialize result
            for (int matrixCol = 0; matrixCol < result.GetNumColumns(); ++matrixCol)
			{
				result(row,matrixCol) = 0.0;
			}
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			// perform multiplication
            for (unsigned int pos = 0; pos < thisColumnVec.size(); ++pos)
			{
				int column = thisColumnVec[pos];
				T value = thisValueVec[pos];
                for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); ++matrixCol)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	return result;
}

//! @brief ... add the scaled other matrix
//! @param rOther ... other matrix
//! @param rFactor ... scalar factor
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::AddScal(const SparseMatrixCSRVector2<T> &rOther, T rFactor)
{
    if (rOther.IsSymmetric())
    {
        AddScal(rOther.AsSparseMatrixCSRVector2Symmetric(), rFactor);
        return;
    }

	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");
	}

    const auto& values = rOther.GetValues();
    const auto& columns = rOther.GetColumns();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int row = 0; row < rOther.GetNumRows(); row++)
    {
        for (unsigned int pos = 0; pos < values[row].size(); pos++)
        {
            this->AddValue(row, columns[row][pos], rFactor*values[row][pos]);
        }
    }
}

//! @brief ... add the scaled other matrix
//! @param rOther ... other matrix
//! @param rFactor ... scalar factor
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::AddScal(const SparseMatrixCSRVector2Symmetric<T> &rOther, T rFactor)
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");
	}
	for (int row = 0; row < rOther.GetNumRows(); row++)
	{
		for (unsigned int pos = 0; pos < rOther.mValues[row].size(); pos++)
		{
			this->AddValue(row, rOther.mColumns[row][pos], rFactor*rOther.mValues[row][pos]);
			if (row!=rOther.mColumns[row][pos])
				this->AddValue(rOther.mColumns[row][pos], row, rFactor*rOther.mValues[row][pos]);
		}
	}
}

//! @brief ... add the scaled other matrix
//! @param rOther ... other matrix
//! @param rFactor ... scalar factor
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::AddScal(const SparseMatrixCSRGeneral<T> &rOther, T rFactor)
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");
	}
	for (int row = 0; row < rOther.GetNumRows(); row++)
	{
		for (int pos = rOther.mRowIndex[row]; pos < rOther.mRowIndex[row + 1]; pos++)
		{
			this->AddValue(row, rOther.mColumns[pos], rFactor*rOther.mValues[pos]);
		}
	}
}

//! @brief ... add the scaled other matrix
//! @param rOther ... other matrix
//! @param rFactor ... scalar factor
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::AddScal(const SparseMatrixCSRSymmetric<T> &rOther, T rFactor)
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");
	}
	for (int row = 0; row < rOther.GetNumRows(); row++)
	{
		for (int pos = rOther.mRowIndex[row]; pos < rOther.mRowIndex[row + 1]; pos++)
		{
			this->AddValue(row, rOther.mColumns[pos], rFactor*rOther.mValues[pos]);
			if (row!=rOther.mColumns[pos])
				this->AddValue(rOther.mColumns[pos], row, rFactor*rOther.mValues[pos]);
		}
	}

}

template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Add_TransA_B_C_Scal(
        const SparseMatrixCSRVector2<T>& rA,
        const SparseMatrixCSRVector2<T>& rB,
        const SparseMatrixCSRVector2<T>& rC, T rScalar)
{
    int resultNumRows = rA.GetNumColumns();
    int resultNumCols = rC.GetNumColumns();

    if (rA.IsSymmetric() or rB.IsSymmetric() or rC.IsSymmetric())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] A, B and C must be in General format.");

    if (rA.GetNumRows() != rB.GetNumRows())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong matrix dimensions for A.T * B");

    if (rC.GetNumRows() != rB.GetNumColumns())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong matrix dimensions for (A.T * B) * C");

    if (resultNumRows != this->GetNumRows() or resultNumCols != this->GetNumColumns())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] matrix dimension mismatch of *this and A.T * B");

    if (this->HasOneBasedIndexing() or rA.HasOneBasedIndexing() or rB.HasOneBasedIndexing())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] all matrices must have zero based indexing.");

    const auto& A = rA.AsSparseMatrixCSRVector2General();
    const auto& B = rB.AsSparseMatrixCSRVector2General();
    const auto& C = rC.AsSparseMatrixCSRVector2General();

    // not easy to parallelize since iCol is not the access row in AddValue
    // loop over columns of transpose of A
    for (int iCol = 0; iCol < A.GetNumRows(); ++iCol)
    {
        const auto& aTValueVec = A.mValues[iCol];
        const auto& aTRowVec = A.mColumns[iCol];

        for (unsigned int iRow = 0; iRow < aTRowVec.size(); ++iRow)
        {
            T aTValue = aTValueVec[iRow];
            int aTRow = aTRowVec[iRow]; // A is transposed


            const auto& bValueVec = B.mValues[iCol];
            const auto& bColumnVec = B.mColumns[iCol];

            for (unsigned int i = 0; i < bColumnVec.size(); ++i)
            {
                T aTbValue = aTValue * bValueVec[i];
                int bCol = bColumnVec[i];
//                tmp.AddValue(aTRow, bCol, aTbValue);

                const auto& cValueVec = C.mValues[bCol];
                const auto& cColumnVec = C.mColumns[bCol];

                for (unsigned int j = 0; j < cColumnVec.size(); ++j)
                {
                    this->AddValue(aTRow, cColumnVec[j], aTbValue*cValueVec[j]*rScalar);
                }
            }
        }
    }
}

template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Sub_TransA_B_Plus_C_D_Scal(
        const SparseMatrixCSRVector2<T>& rA,
        const SparseMatrixCSRVector2<T>& rB,
        const SparseMatrixCSRVector2<T>& rC,
        const SparseMatrixCSRVector2<T>& rD, T rScalar)
{
    int resultNumRows = rA.GetNumColumns();
    int resultNumCols = rB.GetNumColumns();

    if (resultNumRows != rC.GetNumRows() or resultNumCols != rD.GetNumColumns())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] dim(A.T B) != dim (C D)");


    if (resultNumRows != this->GetNumRows() or resultNumCols != this->GetNumColumns())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] matrix dimension mismatch of *this and A.T * B");


    if (this->HasOneBasedIndexing() or rA.HasOneBasedIndexing() or rB.HasOneBasedIndexing() or rC.HasOneBasedIndexing() or rD.HasOneBasedIndexing())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] all matrices must have zero based indexing.");


    /*
     * -= A.T * B
     */

    // not easy to parallelize since iCol is not the access row in AddValue
    // loop over columns of transpose of A
    for (int iCol = 0; iCol < rA.GetNumRows(); ++iCol)
    {
        const auto& aTValueVec = rA.GetValues()[iCol];
        const auto& aTRowVec = rA.GetColumns()[iCol];

        for (unsigned int iRow = 0; iRow < aTRowVec.size(); ++iRow)
        {
            T aTValue = aTValueVec[iRow];
            int aTRow = aTRowVec[iRow]; // A is transposed


            const auto& bValueVec = rB.GetValues()[iCol];
            const auto& bColumnVec = rB.GetColumns()[iCol];

            for (unsigned int i = 0; i < bColumnVec.size(); ++i)
            {
                T bValue = bValueVec[i];
                int bCol = bColumnVec[i];
                this->AddValue(aTRow, bCol, - aTValue * bValue * rScalar);
            }
        }
    }


    /*
     * -= C * D
     */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int cRow = 0; cRow < rC.GetNumRows(); ++cRow)
    {
        const auto& cValueVec = rC.GetValues()[cRow];
        const auto& cColumnVec = rC.GetColumns()[cRow];

        for (unsigned int iCol = 0; iCol < cColumnVec.size(); ++iCol)
        {
            unsigned int cCol = cColumnVec[iCol];
            T cValue = cValueVec[iCol];

            const auto& dColumnVec = rD.GetColumns()[cCol];
            const auto& dValueVec = rD.GetValues()[cCol];


            for (unsigned int i = 0; i < dColumnVec.size(); ++i)
            {
                this->AddValue(cRow, dColumnVec[i], - cValue * dValueVec[i] * rScalar);
            }
        }
    }
}

template<class T>
NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> NuTo::SparseMatrixCSRVector2General<T>::TransMult(const NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix) const
{
	if (this->GetNumRows() != rMatrix.GetNumRows())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> result(this->GetNumColumns(),rMatrix.GetNumColumns());
	if (this->HasOneBasedIndexing())
	{
		// loop over columns of transpose
		for (int column = 0; column < this->GetNumRows(); column++)
		{
			// perform multiplication
			const std::vector<T>& thisValueVec(this->mValues[column]);
			const std::vector<int>& thisColumnVec(this->mColumns[column]);
			for (unsigned int pos = 0; pos < thisValueVec.size(); pos++)
			{
				int row = thisColumnVec[pos] - 1;
				T value = thisValueVec[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	else
	{
		// loop over columns of transpose
		for (int column = 0; column < this->GetNumRows(); column++)
		{
			// perform multiplication
			const std::vector<T>& thisValueVec(this->mValues[column]);
			const std::vector<int>& thisColumnVec(this->mColumns[column]);
			for (unsigned int pos = 0; pos < thisValueVec.size(); pos++)
			{
				int row = thisColumnVec[pos];
				T value = thisValueVec[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	return result;
}

//! @brief ... calculate the transpose of the matrix (transpose row and columns)
//! @return ... transpose of this matrix (sparse csr storage)
template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::Transpose() const
{
	SparseMatrixCSRVector2General<T> transMatrix(this->GetNumColumns(), this->GetNumRows());
	if (this->mOneBasedIndexing)
	{
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
			{
				transMatrix.AddValue(thisColumnVec[pos]-1, row, thisValueVec[pos]);
			}
		}
		transMatrix.SetOneBasedIndexing();
	}
	else
	{
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
			{
				transMatrix.AddValue(thisColumnVec[pos], row, thisValueVec[pos]);
			}
		}
	}
	return transMatrix;
}

template<class T>
T NuTo::SparseMatrixCSRVector2General<T>::Sum() const
{
    T sum = 0;
    for (auto rowValues : this->mValues)
        sum += std::accumulate(rowValues.begin(), rowValues.end(), T(0));
    return sum;
}

template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Gauss(NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
{
    throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] : to be implemented");
}

template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Gauss(SparseMatrixCSRVector2<T>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance)
{
    /*
     * TODO: fix this as a student project or whatever...
     */
    if (rRhs.IsSymmetric())
        throw NuTo::MathException(std::string("[") + __PRETTY_FUNCTION__ + "] : rRhs must be in general format.");

    SparseMatrixCSRGeneral<T> generalThis(this->AsSparseMatrixCSRVector2General());
    SparseMatrixCSRGeneral<T> generalRhs(rRhs.AsSparseMatrixCSRVector2General());

    generalThis.Gauss(generalRhs, rMappingNewToInitialOrdering, rMappingInitialToNewOrdering, rRelativeTolerance);

    (*this) = generalThis;
    rRhs.AsSparseMatrixCSRVector2General() = generalRhs;
}


//! @brief ... reorder columns of the matrix
//! @param rMappingInitialToNewOrdering ... mapping from initial to new ordering
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int iRow = 0; iRow < this->GetNumRows(); ++iRow)
    {
        auto& columns = this->mColumns[iRow];
        auto& values = this->mValues[iRow];

        // renumber columns according to mapping
        for (unsigned int iCol = 0; iCol < columns.size(); ++iCol)
            columns[iCol] = rMappingInitialToNewOrdering[columns[iCol]];

        // sort columns and apply same sorting to the values
        // this "zipped" sorting is not trivially solved by std::sort
        // --> own sorting algorithm.
        // optimized insertion sort, see https://en.wikipedia.org/wiki/Insertion_sort
        for (unsigned i = 1; i < columns.size(); ++i)
        {
            int col = columns[i];
            T val = values[i];

            int j = i - 1;
            while (j >= 0 && columns[j] > col) // aborts if j < 0, so j cannot be unsigned. or change the algorithm...
            {
                columns[j + 1] = columns[j];
                values[j + 1] = values[j];
                --j;
            }
            columns[j + 1] = col;
            values[j + 1] = val;
        }
    }

}

//! @brief ... remove columns from the end of the matrix
//! @param rNumColumn ... number of colums to be removed
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::RemoveLastColumns(unsigned int rNumColumns)
{
    if (this->HasOneBasedIndexing())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] zero based indexing required. Fixing that (if you need it) should be trivial.");

    this->mNumColumns -= rNumColumns;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int iRow = 0; iRow < this->GetNumRows(); ++iRow)
    {
        auto& columns = this->mColumns[iRow];
        auto& values = this->mValues[iRow];

        auto itToFirstColumnToDelete = std::lower_bound(columns.begin(), columns.end(), this->mNumColumns);
        auto itToFirstValueToDelete = values.begin() + std::distance(columns.begin(), itToFirstColumnToDelete);

        columns.erase(itToFirstColumnToDelete, columns.end());
        values.erase(itToFirstValueToDelete, values.end());
    }



}
//! @brief ... Concatenate columns from another matrix
//! @param rOther ... other matrix with same number of rows
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::ConcatenateColumns(const SparseMatrixCSRVector2General<T>& rOther)
{
/*    //! @brief value of nonzero matrix entries
    std::vector<std::vector<T> > mValues;
    //! @brief columns of nonzero matrix entries
    std::vector<std::vector<int> > mColumns;
    //! @brief ... number of columns
    int mNumColumns;
*/
    if (this->mOneBasedIndexing!=rOther.mOneBasedIndexing)
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] index base (0 or 1) should be identical for both matrices");
    if (this->mValues.size()!=rOther.mValues.size())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] number of rows has to be identical for both matrices.");

    for (unsigned int theRow=0; theRow<this->mValues.size(); theRow++)
    {
        this->mValues[theRow].insert(this->mValues[theRow].end(), rOther.mValues[theRow].begin(), rOther.mValues[theRow].end());
        unsigned int oldSize(this->mColumns[theRow].size());
        this->mColumns[theRow].insert(this->mColumns[theRow].end(), rOther.mColumns[theRow].begin(), rOther.mColumns[theRow].end());
        //add to all columns of the newly added entries the columns of the original matrix
        for (unsigned int theCol=oldSize; theCol<this->mColumns[theRow].size(); theCol++)
        {
            this->mColumns[theRow][theCol]+=this->mNumColumns;
        }
    }
    this->mNumColumns+=rOther.mNumColumns;
}

//! @brief ... Concatenate rows from another matrix
//! @param rOther ... other matrix with same number of columns
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::ConcatenateRows(const SparseMatrixCSRVector2General<T>& rOther)
{
    if (this->mOneBasedIndexing!=rOther.mOneBasedIndexing)
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] index base (0 or 1) should be identical for both matrices");
    if (this->mNumColumns!=rOther.mNumColumns)
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] number of columns has to be identical for both matrices.");

    this->mValues.insert(this->mValues.end(), rOther.mValues.begin(), rOther.mValues.end());
    this->mColumns.insert(this->mColumns.end(), rOther.mColumns.begin(), rOther.mColumns.end());
}

//! @brief ... resize matrix
//! @param rNumRows_ ... number of rows
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Resize(int rNumRows, int rNumColumns)
{
    // check for overflow
    assert(rNumRows < INT_MAX);
    assert(rNumRows >= 0);

    //no resize, since the reserved size is only decreased with a copy of a new object
    this->mValues  = std::vector<std::vector<T> >(rNumRows);
    this->mColumns = std::vector<std::vector<int > >(rNumRows);

    mNumColumns = rNumColumns;
}

//! @brief ... print info about the object
template<class T>
void NuTo::SparseMatrixCSRVector2General<T>::Info() const
{
	std::cout << "number of rows: " << this->mValues.size() << std::endl;
	std::cout << "number of columns: " << this->mNumColumns << std::endl;
    for (unsigned int row = 0; row < this->mValues.size(); row++)
    {
        std::cout << "row " << row << ": ";
        for (unsigned int col_count=0; col_count<this->mValues[row].size(); col_count++)
        	if (this->mOneBasedIndexing)
    	        std::cout << "col " << this->mColumns[row][col_count]-1 << " val " << this->mValues[row][col_count] << "    ";
        	else
	            std::cout << "col " << this->mColumns[row][col_count] << " val " << this->mValues[row][col_count] << " ";
        std::cout << std::endl;
    }
}

template<class T>
NuTo::SparseMatrixCSRVector2General<T> NuTo::SparseMatrixCSRVector2General<T>::Random(int rNumRows, int rNumColumns, double rDensity, int rSeed)
{
    SparseMatrixCSRVector2General<T> matrix(rNumRows, rNumColumns);
    Matrix<T>::FillMatrixRandom(matrix, rDensity, rSeed);
    return std::move(matrix);
}

template <class T>
NuTo::SparseMatrixCSRVector2General<T>& NuTo::SparseMatrixCSRVector2General<T>::AsSparseMatrixCSRVector2General()
{
	return *this;
}

template <class T>
const NuTo::SparseMatrixCSRVector2General<T>& NuTo::SparseMatrixCSRVector2General<T>::AsSparseMatrixCSRVector2General()const
{
	return *this;
}

