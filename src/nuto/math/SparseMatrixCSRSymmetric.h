// $Id$

#ifndef SPARSE_MATRIX_CSR_SYMMETRIC_H
#define SPARSE_MATRIX_CSR_SYMMETRIC_H

#include "nuto/math/MathException.h"

#include "nuto/math/SparseMatrixCSRSymmetric_Def.h"

#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

//! @brief ... constructor
//! @param rDimension_ ... dimension (number of rows and number of columns) of square matrix
//! @param rNumReserveEntries_ ... number of entries for which memory is reserved (optional)
template <class T>
NuTo::SparseMatrixCSRSymmetric<T>::SparseMatrixCSRSymmetric(int rDimension, int rNumReserveEntries) : SparseMatrixCSR<T>(rDimension, rNumReserveEntries)
{
}

//! @brief ... returns whether the matrix is symmetric or unsymmetric
//! @return true if the matrix is symmetric and false if the matrix is unsymmetric
template <class T>
bool NuTo::SparseMatrixCSRSymmetric<T>::IsSymmetric() const
{
	return true;
}

//! @brief ... returns the number of columns
//! @return number of columns
template <class T>
int NuTo::SparseMatrixCSRSymmetric<T>::GetNumColumns() const
{
	return this->mRowIndex.size() -  1;
}

//! @brief ... add nonzero entry to matrix
//! @param rRow ... row of the nonzero entry (zero based indexing!!!)
//! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
//! @param rValue ... value of the nonzero entry
template <class T>
void NuTo::SparseMatrixCSRSymmetric<T>::AddEntry(int rRow, int rColumn, T rValue)
{
	// check for overflow
	assert(rRow < INT_MAX);
	assert(rColumn < INT_MAX);

	// check bounds
	if (rRow >= (int)this->mRowIndex.size() - 1)
	{
		throw MathException("[SparseMatrixCSRSymmetric::addEntry] row index is out of bounds.");
	}
	if (rColumn >= (int)this->mRowIndex.size() - 1)
	{
		throw MathException("[SparseMatrixCSRSymmetric::addEntry] column index is out of bounds.");
	}
	if (rColumn < rRow)
	{
		throw MathException("[SparseMatrixCSRSymmetric::addEntry] upper triangle is stored for symmetric matrices.");
	}
	if (this->mOneBasedIndexing)
	{
		rColumn++;

		// find position in matrix
		int pos = this->mRowIndex[rRow] - 1;
		for (; (pos < this->mRowIndex[rRow + 1] - 1) && this->mColumns[pos] < static_cast<int>(rColumn); pos++);

		// add value
		if ((pos == this->mRowIndex[rRow + 1] - 1) || (static_cast<int>(rColumn) != this->mColumns[pos]))
		{
			// insert new value
			this->mColumns.insert(this->mColumns.begin() + pos, rColumn);
			this->mValues.insert(this->mValues.begin() + pos, rValue);
			for (unsigned int row_count = rRow + 1; row_count < this->mRowIndex.size(); row_count++)
			{
				this->mRowIndex[row_count] += 1;
			}
		}
		else
		{
			// add to existing value
			this->mValues[pos] += rValue;
		}
	}
	else // zero based indexing
	{
		// find position in matrix
		int pos = this->mRowIndex[rRow];
		for (; (pos < this->mRowIndex[rRow + 1]) && this->mColumns[pos] < static_cast<int>(rColumn); pos++);

		// add value
		if ((pos == this->mRowIndex[rRow + 1]) || (static_cast<int>(rColumn) != this->mColumns[pos]))
		{
			// insert new value
			this->mColumns.insert(this->mColumns.begin() + pos, rColumn);
			this->mValues.insert(this->mValues.begin() + pos, rValue);
			for (unsigned int row_count = rRow + 1; row_count < this->mRowIndex.size(); row_count++)
			{
				this->mRowIndex[row_count] += 1;
			}
		}
		else
		{
			// add to existing value
			this->mValues[pos] += rValue;
		}
	}
}

//! @brief ... print info about the object
template <class T>
void NuTo::SparseMatrixCSRSymmetric<T>::Info() const
{
	SparseMatrixCSR<T>::Info();
}

//! @brief ... write nonzero matrix entries into a full matrix
//! @param rFullMatrix ... the full matrix
template <class T>
void NuTo::SparseMatrixCSRSymmetric<T>::WriteEntriesToFullMatrix(FullMatrix<T>& rFullMatrix) const
{
	std::vector<int>::const_iterator columnIterator = this->mColumns.begin();
	typename std::vector<T>::const_iterator valueIterator = this->mValues.begin();
	if (this->mOneBasedIndexing)
	{
		unsigned int row = 0;
		while (row < this->mRowIndex.size() - 1)
		{
			for (int entry_count = this->mRowIndex[row]; entry_count < this->mRowIndex[row+1]; entry_count++)
			{
				rFullMatrix(row, (*columnIterator) - 1) = *valueIterator;
				rFullMatrix((*columnIterator) - 1, row) = *valueIterator;
				columnIterator++;
				valueIterator++;
			}
			row++;
		}
	}
	else
	{
		unsigned int row = 0;
		while (row < this->mRowIndex.size() - 1)
		{
			for (int entry_count = this->mRowIndex[row]; entry_count < this->mRowIndex[row+1]; entry_count++)
			{
				// set upper triangle and diagonal entries
				rFullMatrix(row, *columnIterator) = *valueIterator;
				// set lower triangle entries
				if (static_cast<int>(row) != *columnIterator)
				{
					rFullMatrix(*columnIterator, row) = *valueIterator;
				}
				columnIterator++;
				valueIterator++;
			}
			row++;
		}
	}
}

#endif // SPARSE_MATRIX_CSR_SYMMETRIC_H
