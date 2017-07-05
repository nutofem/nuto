
#pragma once

#include "base/Exception.h"

#include "math/SparseMatrixCSRSymmetric_Def.h"

#include "math/SparseMatrixCSRVector2Symmetric.h"

//! @brief ... constructor
//! @param rDimension_ ... dimension (number of rows and number of columns) of square matrix
//! @param rNumReserveEntries_ ... number of entries for which memory is reserved (optional)
template <class T>
NuTo::SparseMatrixCSRSymmetric<T>::SparseMatrixCSRSymmetric(int rDimension, int rNumReserveEntries)
    : SparseMatrixCSR<T>(rDimension, rNumReserveEntries)
{
}

//! @brief ... constructor
//! @param rDimension_ ... dimension (number of rows and number of columns) of square matrix
//! @param rNumReserveEntries_ ... number of entries for which memory is reserved (optional)
template <class T>
NuTo::SparseMatrixCSRSymmetric<T>::SparseMatrixCSRSymmetric(const SparseMatrixCSRSymmetric<T>& rOther)
    : SparseMatrixCSR<T>(rOther)
{
}

//! @brief ... constructor
//! @param rDimension_ ... dimension (number of rows and number of columns) of square matrix
//! @param rNumReserveEntries_ ... number of entries for which memory is reserved (optional)
template <class T>
NuTo::SparseMatrixCSRSymmetric<T>::SparseMatrixCSRSymmetric(const SparseMatrixCSRVector2Symmetric<T>& rOther)
    : SparseMatrixCSR<T>(rOther.GetNumRows())
{
    int numEntries = rOther.GetNumEntries();
    this->mValues.resize(numEntries);
    this->mColumns.resize(numEntries);

    int theEntry(0);
    assert(rOther.mValues.size() == rOther.mColumns.size());
    for (int row = 0; row < rOther.GetNumRows(); row++)
    {
        this->mRowIndex[row] = theEntry + this->mOneBasedIndexing;
        for (unsigned int colCount = 0; colCount < rOther.mValues[row].size(); colCount++)
        {
            assert(rOther.mValues[row].size() == rOther.mColumns[row].size());
            assert(theEntry < numEntries);
            this->mValues[theEntry] = rOther.mValues[row][colCount];
            this->mColumns[theEntry] =
                    rOther.mColumns[row][colCount] - rOther.mOneBasedIndexing + this->mOneBasedIndexing;
            theEntry++;
        }
    }
    this->mRowIndex[rOther.GetNumRows()] = theEntry + this->mOneBasedIndexing;
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
    return this->mRowIndex.size() - 1;
}

//! @brief ... add nonzero entry to matrix
//! @param rRow ... row of the nonzero entry (zero based indexing!!!)
//! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
//! @param rValue ... value of the nonzero entry
template <class T>
void NuTo::SparseMatrixCSRSymmetric<T>::AddValue(int rRow, int rColumn, const T& rValue)
{
	// check for overflow
	assert(rRow < INT_MAX);
	assert(rColumn < INT_MAX);

	// check bounds
	if (rRow >= (int)this->mRowIndex.size() - 1)
	{
		throw Exception("[SparseMatrixCSRSymmetric::addEntry] row index is out of bounds.");
	}
	if (rColumn >= (int)this->mRowIndex.size() - 1)
	{
		throw Exception("[SparseMatrixCSRSymmetric::addEntry] column index is out of bounds.");
	}
	if (rColumn < rRow)
	{
		throw Exception("[SparseMatrixCSRSymmetric::addEntry] upper triangle is stored for symmetric matrices.");
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

//! @brief ... return the matrix type
template <class T>
NuTo::eSparseMatrixType NuTo::SparseMatrixCSRSymmetric<T>::GetSparseMatrixType() const
{
    return NuTo::eSparseMatrixType::CSRSYMMETRIC;
}

//! @brief ... print info about the object
template <class T>
void NuTo::SparseMatrixCSRSymmetric<T>::Info() const
{
    SparseMatrixCSR<T>::Info();
}

//! @brief ... write nonzero matrix entries into a matrix
//! @param rFullMatrix ... the full matrix
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> NuTo::SparseMatrixCSRSymmetric<T>::ConvertToFullMatrix() const
{
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> m =
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(this->GetNumRows(), this->GetNumColumns());
    std::vector<int>::const_iterator columnIterator = this->mColumns.begin();
    typename std::vector<T>::const_iterator valueIterator = this->mValues.begin();
    if (this->mOneBasedIndexing)
    {
        unsigned int row = 0;
        while (row < this->mRowIndex.size() - 1)
        {
            for (int entry_count = this->mRowIndex[row]; entry_count < this->mRowIndex[row + 1]; entry_count++)
            {
                m(row, (*columnIterator) - 1) = *valueIterator;
                if (static_cast<int>(row) != (*columnIterator) - 1)
                    m((*columnIterator) - 1) = row, *valueIterator;
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
            for (int entry_count = this->mRowIndex[row]; entry_count < this->mRowIndex[row + 1]; entry_count++)
            {
                // set upper triangle and diagonal entries
                m(row, *columnIterator) = *valueIterator;
                // set lower triangle entries
                if (static_cast<int>(row) != *columnIterator)
                {
                    m(*columnIterator, row) = *valueIterator;
                }
                columnIterator++;
                valueIterator++;
            }
            row++;
        }
    }
    return m;
}

//! @brief ... resize matrix
//! @param rNumRows_ ... number of rows
template <class T>
void NuTo::SparseMatrixCSRSymmetric<T>::Resize(int rNumRows)
{
    // check for overflow
    assert(rNumRows < INT_MAX);
    assert(rNumRows >= 0);

    // resize
    SparseMatrixCSR<T>::Resize(rNumRows);
}

//! @brief ... resize matrix
//! @param rNumRows_ ... number of rows
//! @param rNumColumns_ ... number of columns
template <class T>
void NuTo::SparseMatrixCSRSymmetric<T>::Resize(int rNumRows_, int rNumColumns_)
{
	// check for overflow
	assert(rNumColumns_ < INT_MAX);
	assert(rNumColumns_ >= 0);
	if (rNumRows_!=rNumColumns_)
		throw Exception("[SparseMatrixCSRSymmetric::Resize] number of rows and column has to be identical for symmetric matrices.");

	// resize
	SparseMatrixCSR<T>::Resize(rNumRows_);
}

template <class T>
NuTo::SparseMatrixCSRSymmetric<T>& NuTo::SparseMatrixCSRSymmetric<T>::AsSparseMatrixCSRSymmetric()
{
    return *this;
}

template <class T>
const NuTo::SparseMatrixCSRSymmetric<T>& NuTo::SparseMatrixCSRSymmetric<T>::AsSparseMatrixCSRSymmetric() const
{
    return *this;
}
