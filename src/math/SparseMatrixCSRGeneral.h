// $Id$

#pragma once
#include "math/SparseMatrixCSRGeneral_Def.h"

#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseMatrixEnum.h"
#include "math/MathException.h"

//! @brief ... constructor
//! @param rNumRows_ ... number of rows
//! @param rNumColumns_ ... number of columns
//! @param rNumReserveEntries_ ... number of entries for which memory is reserved (optional)
template <class T>
NuTo::SparseMatrixCSRGeneral<T>::SparseMatrixCSRGeneral(int rNumRows_, int rNumColumns_, unsigned int rNumReserveEntries_) : NuTo::SparseMatrixCSR<T>(rNumRows_, rNumReserveEntries_)
{
	// check for overflow
	assert(rNumColumns_ < INT_MAX);
	assert(rNumColumns_ >= 0);
	this->mNumColumns = rNumColumns_;
}

//! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a predefined tolerance)
//! @param rFullMatrix ... input matrix (full storage)
//! @param rAbsoluteTolerance ... absolute tolerance
//! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance * max(abs(rMatrixEntry))
template <class T>
NuTo::SparseMatrixCSRGeneral<T>::SparseMatrixCSRGeneral(const NuTo::SparseMatrixCSRVector2General<T>& rCSR2Matrix) : NuTo::SparseMatrixCSR<T>(rCSR2Matrix.GetNumRows(),rCSR2Matrix.GetNumColumns())
{
	//allocate space
	int numEntries = rCSR2Matrix.GetNumEntries();
	this->mValues.resize(numEntries);
	this->mColumns.resize(numEntries);
	this->mRowIndex.resize(rCSR2Matrix.GetNumRows()+1);
	int globalPos(0);
	this->mOneBasedIndexing = rCSR2Matrix.HasOneBasedIndexing();
	for (int row=0; row<rCSR2Matrix.GetNumRows(); row++)
	{
		this->mRowIndex[row] = globalPos+this->mOneBasedIndexing;
		const std::vector<T>& thisValueVec(rCSR2Matrix.mValues[row]);
		const std::vector<int>& thisColumnVec(rCSR2Matrix.mColumns[row]);
		for (unsigned int the_entry=0; the_entry < thisValueVec.size(); the_entry++, globalPos++)
		{
			this->mValues[globalPos]=thisValueVec[the_entry];
			this->mColumns[globalPos]=thisColumnVec[the_entry];
		}
	}
	this->mRowIndex[rCSR2Matrix.GetNumRows()] = globalPos+this->mOneBasedIndexing;
	mNumColumns = rCSR2Matrix.GetNumColumns();
}

//! @brief ... resize matrix
//! @param rNumRows_ ... number of rows
//! @param rNumColumns_ ... number of columns
template <class T>
void NuTo::SparseMatrixCSRGeneral<T>::Resize(int rNumRows_, int rNumColumns_)
{
	// check for overflow
	assert(rNumColumns_ < INT_MAX);
	assert(rNumColumns_ >= 0);

	// resize
	SparseMatrixCSR<T>::Resize(rNumRows_);
	this->mNumColumns = rNumColumns_;
}

//! @brief ... returns whether the matrix is symmetric or unsymmetric
//! @return true if the matrix is symmetric and false if the matrix is unsymmetric
template <class T>
bool NuTo::SparseMatrixCSRGeneral<T>::IsSymmetric() const
{
	return false;
}

//! @brief ... returns the number of columns
//! @return number of columns
template <class T>
int NuTo::SparseMatrixCSRGeneral<T>::GetNumColumns() const
{
	return this->mNumColumns;
}

//! @brief ... add nonzero entry to matrix
//! @param rRow ... row of the nonzero entry (zero based indexing!!!)
//! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
//! @param rValue ... value of the nonzero entry
template <class T>
void NuTo::SparseMatrixCSRGeneral<T>::AddValue(int rRow, int rColumn, const T& rValue)
{
	// check for overflow
    assert(rRow < INT_MAX);
    assert(rColumn < INT_MAX);
    assert(rRow >= 0);
    assert(rColumn >= 0);

    // check bounds
    if (rRow >= (int)this->mRowIndex.size() - 1 || rRow<0)
    {
        throw MathException("[SparseMatrixCSRGeneral::addEntry] row index is out of bounds.");
    }
    if (rColumn >= this->mNumColumns || rColumn<0)
    {
        throw MathException("[SparseMatrixCSRGeneral::addEntry] column index is out of bounds.");
    }

    if (this->mOneBasedIndexing)
    {
        rColumn++;

        // find position in matrix
        int pos = this->mRowIndex[rRow] - 1;
        for (; (pos < this->mRowIndex[rRow + 1] - 1) && this->mColumns[pos] < static_cast<int>(rColumn); pos++)
            ;

        // add value
        if ((pos == this->mRowIndex[rRow + 1] - 1) || (static_cast<int>(rColumn) != this->mColumns[pos])
           )
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
        for (; (pos < this->mRowIndex[rRow + 1]) && this->mColumns[pos] < static_cast<int>(rColumn); pos++)
            ;

        // add value
        if ((pos == this->mRowIndex[rRow + 1]) || (static_cast<int>(rColumn) != this->mColumns[pos])
           )
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
template<class T>
NuTo::eSparseMatrixType NuTo::SparseMatrixCSRGeneral<T>::GetSparseMatrixType()const
{
    return NuTo::eSparseMatrixType::CSRGENERAL;
}

//! @brief ... print info about the object
template <class T>
void NuTo::SparseMatrixCSRGeneral<T>::Info() const
{
	std::cout << "number of columns: " << this->mNumColumns << std::endl;
	SparseMatrixCSR<T>::Info();
}

template <class T>
Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> NuTo::SparseMatrixCSRGeneral<T>::ConvertToFullMatrix() const
{
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> m = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(this->GetNumRows(), this->GetNumColumns());
	std::vector<int>::const_iterator columnIterator = this->mColumns.begin();
	typename std::vector<T>::const_iterator valueIterator = this->mValues.begin();
	if (this->mOneBasedIndexing)
	{
		unsigned int row = 0;
		while (row < this->mRowIndex.size() - 1)
		{
			for (int entry_count = this->mRowIndex[row]; entry_count < this->mRowIndex[row+1]; entry_count++)
			{
				m(row, (*columnIterator) - 1) = *valueIterator;
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
				m(row, *columnIterator) = *valueIterator;
				columnIterator++;
				valueIterator++;
			}
			row++;
		}
	}
    return m;
}


//! @brief ... add two matrices
//! @param rOther ... general sparse matrix stored in the CSR format
//! @return general sparse matrix stored in the CSR format
template <class T>
NuTo::SparseMatrixCSRGeneral<T> NuTo::SparseMatrixCSRGeneral<T>::operator+ (const NuTo::SparseMatrixCSRGeneral<T> &rOther ) const
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRGeneral::operator+] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRGeneral::operator+] both matrices must have zero based indexing.");
	}
	SparseMatrixCSRGeneral<T> result(this->GetNumRows(), this->GetNumColumns());
	for (int row = 0; row < this->GetNumRows(); row++)
	{
		int thisPos = this->mRowIndex[row];
		int otherPos = rOther.mRowIndex[row];
		while ((thisPos < this->mRowIndex[row + 1]) || (otherPos < rOther.mRowIndex[row + 1]))
		{
			int thisColumn;
			if (thisPos < this->mRowIndex[row + 1])
			{
				thisColumn = this->mColumns[thisPos];
			}
			else
			{
				// no additional entries in this row (set column to an invalid value)
				thisColumn = this->GetNumColumns();
			}

			int otherColumn;
			if (otherPos < rOther.mRowIndex[row + 1])
			{
				otherColumn = rOther.mColumns[otherPos];
			}
			else
			{
				// no additional entries in this row (set column to an invalid value)
				otherColumn = rOther.GetNumColumns();
			}

			int resultColumn;
			T resultValue;
			if (thisColumn < otherColumn)
			{
				resultColumn = thisColumn;
				resultValue = this->mValues[thisPos];
				thisPos++;
			}
			else if (otherColumn < thisColumn)
			{
				resultColumn = otherColumn;
				resultValue = rOther.mValues[otherPos];
				otherPos++;
			}
			else
			{
				resultColumn = thisColumn;
				resultValue = this->mValues[thisPos] + rOther.mValues[otherPos];
				thisPos++;
				otherPos++;
			}
			result.mColumns.push_back(resultColumn);
			result.mValues.push_back(resultValue);
		}
		assert(result.mColumns.size() == result.mValues.size());
		result.mRowIndex[row + 1] = result.mColumns.size();
	}
	return result;
}

//! @brief ... subtract two matrices
//! @param rOther ... general sparse matrix stored in the CSR format
//! @return general sparse matrix stored in the CSR format
template <class T>
NuTo::SparseMatrixCSRGeneral<T> NuTo::SparseMatrixCSRGeneral<T>::operator- (const NuTo::SparseMatrixCSRGeneral<T> &rOther ) const
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRGeneral::operator+] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRGeneral::operator+] both matrices must have zero based indexing.");
	}
	SparseMatrixCSRGeneral<T> result(this->GetNumRows(), this->GetNumColumns());
	for (int row = 0; row < this->GetNumRows(); row++)
	{
		int thisPos = this->mRowIndex[row];
		int otherPos = rOther.mRowIndex[row];
		while ((thisPos < this->mRowIndex[row + 1]) || (otherPos < rOther.mRowIndex[row + 1]))
		{
			int thisColumn;
			if (thisPos < this->mRowIndex[row + 1])
			{
				thisColumn = this->mColumns[thisPos];
			}
			else
			{
				// no additional entries in this row (set column to an invalid value)
				thisColumn = this->GetNumColumns();
			}

			int otherColumn;
			if (otherPos < rOther.mRowIndex[row + 1])
			{
				otherColumn = rOther.mColumns[otherPos];
			}
			else
			{
				// no additional entries in this row (set column to an invalid value)
				otherColumn = rOther.GetNumColumns();
			}

			int resultColumn;
			T resultValue;
			if (thisColumn < otherColumn)
			{
				resultColumn = thisColumn;
				resultValue = this->mValues[thisPos];
				thisPos++;
			}
			else if (otherColumn < thisColumn)
			{
				resultColumn = otherColumn;
				resultValue = -rOther.mValues[otherPos];
				otherPos++;
			}
			else
			{
				resultColumn = thisColumn;
				resultValue = this->mValues[thisPos] - rOther.mValues[otherPos];
				thisPos++;
				otherPos++;
			}
			result.mColumns.push_back(resultColumn);
			result.mValues.push_back(-resultValue);
		}
		assert(result.mColumns.size() == result.mValues.size());
		result.mRowIndex[row + 1] = result.mColumns.size();
	}
	return result;
}

//! @brief ... subtract two matrices
//! @param rOther ... general sparse matrix stored in the CSR format
//! @return reference to this matrix
template <class T>
NuTo::SparseMatrixCSRGeneral<T>& NuTo::SparseMatrixCSRGeneral<T>::operator-=  (const NuTo::SparseMatrixCSRGeneral<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRGeneral::operator-=] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRGeneral::operator-=] both matrices must have zero based indexing.");
	}
	for (int row = 0; row < rOther.GetNumRows(); row++)
	{
		for (int pos = rOther.mRowIndex[row]; pos < rOther.mRowIndex[row + 1]; pos++)
		{
			this->AddValue(row, rOther.mColumns[pos], - rOther.mValues[pos]);
		}
	}
	return *this;
}


//! @brief ... add two matrices
//! @param rOther ... general sparse matrix stored in the CSR format
//! @return reference to this matrix
template <class T>
NuTo::SparseMatrixCSRGeneral<T>& NuTo::SparseMatrixCSRGeneral<T>::operator+=  (const NuTo::SparseMatrixCSRGeneral<T> &rOther )
{
	if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
	{
		throw MathException("[SparseMatrixCSRGeneral::operator+=] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRGeneral::operator+=] both matrices must have zero based indexing.");
	}
	for (int row = 0; row < rOther.GetNumRows(); row++)
	{
		for (int pos = rOther.mRowIndex[row]; pos < rOther.mRowIndex[row + 1]; pos++)
		{
			this->AddValue(row, rOther.mColumns[pos], rOther.mValues[pos]);
		}
	}
	return *this;
}

//! @brief ... matrix - matrix multiplication
//! @param rOther ... general sparse matrix stored in the CSR format
//! @return general sparse matrix stored in the CSR format
template <class T>
NuTo::SparseMatrixCSRGeneral<T> NuTo::SparseMatrixCSRGeneral<T>::operator* (const NuTo::SparseMatrixCSRGeneral<T> &rOther ) const
{
	if (this->GetNumColumns() != rOther.GetNumRows())
	{
		throw MathException("[SparseMatrixCSRGeneral::operator*] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException("[SparseMatrixCSRGeneral::operator*] both matrices must have zero based indexing.");
	}
	SparseMatrixCSRGeneral<T> result(this->GetNumRows(), rOther.GetNumColumns());

	for (int thisRow = 0; thisRow < this->GetNumRows(); thisRow++)
	{
		for (int thisPos = this->mRowIndex[thisRow]; thisPos < this->mRowIndex[thisRow + 1]; thisPos++)
		{
			unsigned int thisColumn = this->mColumns[thisPos];
			T thisValue = this->mValues[thisPos];

			for (int otherPos = rOther.mRowIndex[thisColumn]; otherPos < rOther.mRowIndex[thisColumn + 1]; otherPos++)
			{
				unsigned int otherColumn = rOther.mColumns[otherPos];
				T otherValue = rOther.mValues[otherPos];
				result.AddValue(thisRow, otherColumn, thisValue * otherValue);
			}
		}
	}
	return result;
}

//! @brief ... multiplies the matrix with an scalar value
//! @param rOther ... scalar value
//! @return ... the multiplied matrix (sparse csr storage)
template <class T>
NuTo::SparseMatrixCSRGeneral<T> NuTo::SparseMatrixCSRGeneral<T>::operator* (const T &rOther ) const
{
	SparseMatrixCSRGeneral<T> result(*this);
	BOOST_FOREACH( T &val, result.mValues )
		val *= rOther;
	return result;
}

//! @brief ... multiply sparse matrix with a full matrix
//! @param rFullMatrix ... full matrix which is multiplied with the sparse matrix
//! @return ... full matrix
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> NuTo::SparseMatrixCSRGeneral<T>::operator* (const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &rMatrix) const
{

	if (this->GetNumColumns() != rMatrix.rows())
	{
		throw MathException("[SparseMatrixCSRGeneral::operator*] invalid matrix dimensions.");
	}
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(this->GetNumRows(),rMatrix.cols());
	if (this->HasOneBasedIndexing())
	{
		// loop over rows
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			// initialize result
			for (int matrixCol = 0; matrixCol < result.cols(); matrixCol++)
			{
				result(row,matrixCol) = 0.0;
			}
			// perform multiplication
			for (int pos = this->mRowIndex[row] - 1; pos < this->mRowIndex[row + 1] - 1; pos++)
			{
				int column = this->mColumns[pos] - 1;
				T value = this->mValues[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.cols(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	else
	{
		// loop over rows
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			// initialize result
			for (int matrixCol = 0; matrixCol < result.cols(); matrixCol++)
			{
				result(row,matrixCol) = 0.0;
			}
			// perform multiplication
			for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
			{
				int column = this->mColumns[pos];
				T value = this->mValues[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.cols(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
			}
		}
	}
	return result;
}

template <class T>
Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> NuTo::SparseMatrixCSRGeneral<T>::TransMult(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix) const
{
	if (this->GetNumRows() != rMatrix.rows())
	{
		throw MathException(__PRETTY_FUNCTION__, "Invalid matrix dimensions.");
	}

	Eigen::Matrix<T, Eigen::Dynamic,Eigen::Dynamic> result(this->GetNumColumns(),rMatrix.cols());
	if (this->HasOneBasedIndexing())
	{
		// loop over columns of transpose
		for (int column = 0; column < this->GetNumRows(); column++)
		{
			// perform multiplication
			for (int pos = this->mRowIndex[column] - 1; pos < this->mRowIndex[column + 1] - 1; pos++)
			{
				int row = this->mColumns[pos] - 1;
				T value = this->mValues[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.cols(); matrixCol++)
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
			for (int pos = this->mRowIndex[column]; pos < this->mRowIndex[column + 1]; pos++)
			{
				int row = this->mColumns[pos];
				T value = this->mValues[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.cols(); matrixCol++)
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
template <class T>
NuTo::SparseMatrixCSRGeneral<T> NuTo::SparseMatrixCSRGeneral<T>::Transpose() const
{
	SparseMatrixCSRGeneral<T> transMatrix(this->GetNumColumns(), this->GetNumRows());
	if (this->mOneBasedIndexing)
	{
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			for (int pos = this->mRowIndex[row] - 1; pos < this->mRowIndex[row + 1] - 1; pos++)
			{
				transMatrix.AddValue(this->mColumns[pos] - 1, row, this->mValues[pos]);
			}
		}
		transMatrix.SetOneBasedIndexing();
	}
	else
	{
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
			{
				transMatrix.AddValue(this->mColumns[pos], row, this->mValues[pos]);
			}
		}
	}
	return transMatrix;
}

//! @brief ... remove entry from matrix
//! @param rRow ... row index (zero based indexing!!!)
//! @param rColumn ... column index (zero based indexing!!!)
template <class T>
void NuTo::SparseMatrixCSRGeneral<T>::RemoveEntry(int rRow, int rColumn)
{
	// check bounds
	if (rRow >= (int)this->mRowIndex.size() - 1 || rRow<0)
	{
		throw MathException("[SparseMatrixCSRGeneral::RemoveEntry] row index is out of bounds.");
	}
	if (rColumn >= this->mNumColumns || rColumn<0)
	{
		throw MathException("[SparseMatrixCSRGeneral::RemoveEntry] column index is out of bounds.");
	}
	if (this->mOneBasedIndexing)
	{
		rColumn++;
		for (int pos = this->mRowIndex[rRow] - 1; pos < this->mRowIndex[rRow + 1] - 1; pos++)
		{
			if (this->mColumns[pos] > rColumn)
			{
				break;
			}
			else if (this->mColumns[pos] == rColumn)
			{
				this->mColumns.erase(this->mColumns.begin() + pos);
				this->mValues.erase(this->mValues.begin() + pos);
				for (unsigned int rowCount = rRow + 1; rowCount < this->mRowIndex.size(); rowCount++)
				{
					this->mRowIndex[rowCount] -= 1;
				}
				break;
			}
		}
	}
	else
	{
		for (int pos = this->mRowIndex[rRow]; pos < this->mRowIndex[rRow + 1]; pos++)
		{
			if (this->mColumns[pos] > rColumn)
			{
				break;
			}
			else if (this->mColumns[pos] == rColumn)
			{
				this->mColumns.erase(this->mColumns.begin() + pos);
				this->mValues.erase(this->mValues.begin() + pos);
				for (unsigned int rowCount = rRow + 1; rowCount < this->mRowIndex.size(); rowCount++)
				{
					this->mRowIndex[rowCount] -= 1;
				}
				break;
			}
		}
	}
}

//! @brief ... remove columns from the end of the matrix
//! @param rNumColumn ... number of colums to be removed
template <class T>
void NuTo::SparseMatrixCSRGeneral<T>::RemoveLastColumns(unsigned int rNumColumns)
{
	if ((unsigned int)this->mNumColumns > rNumColumns)
	{
		if (this->mOneBasedIndexing)
		{
			for (int row = 0; row < this->GetNumRows(); row++)
			{
				for (int pos = this->mRowIndex[row] - 1; pos < this->mRowIndex[row + 1] - 1; pos++)
				{
					if ((unsigned int)this->mColumns[pos] > this->mNumColumns-rNumColumns)
					{
						this->RemoveEntry(row, this->mColumns[pos] - 1);
					}
				}
			}
		}
		else
		{
			for (int row = 0; row < this->GetNumRows(); row++)
			{
				for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
				{
					if ((unsigned int)this->mColumns[pos] >=  this->mNumColumns-rNumColumns)
					{
						this->RemoveEntry(row, this->mColumns[pos]);
					}
				}
			}
		}
		this->mNumColumns -= rNumColumns;
	}
	else
	{
		this->Resize(0,0);
	}
}

//! @brief ... reorder columns of the matrix
//! @param rMappingInitialToNewOrdering ... mapping fron initial to new ordering
template <class T>
void NuTo::SparseMatrixCSRGeneral<T>::ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering)
{
	for (int row = 0; row < this->GetNumRows(); row++)
	{
		for (int pos = this->mRowIndex[row]; pos < this->mRowIndex[row + 1]; pos++)
		{
			this->mColumns[pos] = rMappingInitialToNewOrdering[this->mColumns[pos]];
		}

		// sort columns (simple bubble sort algorithm)
		int start = this->mRowIndex[row];
		int end = this->mRowIndex[row + 1] - 1;
		bool swapFlag;
		do
		{
			swapFlag = false;
			for (int pos = start; pos < end; pos++)
			{
				if (this->mColumns[pos] > this->mColumns[pos + 1])
				{
					swapFlag=true;
					int tmpInt = this->mColumns[pos];
					this->mColumns[pos] = this->mColumns[pos + 1];
					this->mColumns[pos + 1] = tmpInt;

					T tmpDouble = this->mValues[pos];
					this->mValues[pos] = this->mValues[pos + 1];
					this->mValues[pos + 1] = tmpDouble;
				}
			}
			end--;
		}
		while (swapFlag);
	}

}

template <class T>
NuTo::SparseMatrixCSRGeneral<T>& NuTo::SparseMatrixCSRGeneral<T>::AsSparseMatrixCSRGeneral()
{
	return *this;
}

template <class T>
const NuTo::SparseMatrixCSRGeneral<T>& NuTo::SparseMatrixCSRGeneral<T>::AsSparseMatrixCSRGeneral()const
{
	return *this;
}

