// $Id: SparseMatrixCSRVector2Symmetric.h 235 2010-04-22 09:25:38Z arnold2 $

#pragma once

#include <algorithm>

#include "nuto/math/SparseMatrixEnum.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric_Def.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"


#include "nuto/math/FullMatrix.h"
//#include "nuto/math/SparseMatrixCSRSymmetric_Def.h"
//#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/MathException.h"

//! @brief ... constructor
//! @param rNumRows_ ... number of rows
//! @param rNumColumns_ ... number of columns
template<class T>
NuTo::SparseMatrixCSRVector2Symmetric<T>::SparseMatrixCSRVector2Symmetric(int rNumRows_, int rNumColumns_) :
                             NuTo::SparseMatrixCSRVector2<T>(rNumRows_)
{
    if (rNumRows_!=rNumColumns_)
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] Symmetric matrix must have same number of rows and columns.");
}

//! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a predefined tolerance)
//! @param rFullMatrix ... input matrix (full storage)
//! @param rAbsoluteTolerance ... absolute tolerance
//! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance * max(abs(rMatrixEntry))
template<class T>
NuTo::SparseMatrixCSRVector2Symmetric<T>::SparseMatrixCSRVector2Symmetric(const NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rFullMatrix, double rAbsoluteTolerance, double rRelativeTolerance):
                             NuTo::SparseMatrixCSRVector2<T>(rFullMatrix.GetNumRows())
{
	if (rFullMatrix.GetNumColumns()!=rFullMatrix.GetNumRows())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] Symmetric matrix must have same number of rows and columns.");

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
		for (int col = row; col < rFullMatrix.GetNumColumns(); col++)
		{
			if (std::abs(rFullMatrix(row,col)) > tolerance)
			{
				this->AddValue(row,col,rFullMatrix(row,col));
			}
			if (col!=row && std::abs(rFullMatrix(col,row)-rFullMatrix(row,col))>tolerance)
			{
			    throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] Full matrix is not symmetric.");
			}
		}
	}
}

//! @brief ... create sparse matrix with vector of vector from standard CSR format
//! @param rCSRMatrix ... input matrix (full storage)
template<class T>
NuTo::SparseMatrixCSRVector2Symmetric<T>::SparseMatrixCSRVector2Symmetric(const NuTo::SparseMatrixCSRSymmetric<T>& rCSRMatrix) :
              NuTo::SparseMatrixCSRVector2<T>::SparseMatrixCSRVector2(rCSRMatrix.GetNumRows())
{
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

//! @brief ... returns whether the matrix is symmetric or unsymmetric
//! @return true if the matrix is symmetric and false if the matrix is unsymmetric
template<class T>
bool NuTo::SparseMatrixCSRVector2Symmetric<T>::IsSymmetric() const
{
	return true;
}

//! @brief ... returns the number of columns
//! @return number of columns
template <class T>
int NuTo::SparseMatrixCSRVector2Symmetric<T>::GetNumColumns() const
{
	return (int)this->mValues.size();
}

//! @brief ... add nonzero entry to matrix
//! @param rRow ... row of the nonzero entry (zero based indexing!!!)
//! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
//! @param rValue ... value of the nonzero entry
template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::AddValue(int rRow, int rColumn, const T& rValue)
{
    // check bounds via asserts
    assert(rRow < (int)this->mValues.size()     && "row index is out of bounds.");
    assert(rRow >= 0                            && "row index is out of bounds.");

    assert(rColumn < (int)this->mValues.size()  && "column index is out of bounds.");
    assert(rColumn >= 0                         && "column index is out of bounds.");

    assert(rColumn >= rRow                      && "(row,column) not in upper triangle.");

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
NuTo::eSparseMatrixType NuTo::SparseMatrixCSRVector2Symmetric<T>::GetSparseMatrixType()const
{
    return NuTo::eSparseMatrixType::CSRVECTOR2SYMMETRIC;
}

//! @brief ... import matrix from slang object stored in  a text file
//! @param rFileName ... file name
template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::ImportFromSLangText(const char* rFileName)
{
	throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] to be implemented.");
}

//! @brief ... write nonzero matrix entries into a matrix
//! @param rMatrix ... the matrix
template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::WriteEntriesToMatrix(Matrix<T>& rMatrix) const
{
	rMatrix.Resize(this->GetNumRows(), this->GetNumColumns());
	if (this->mOneBasedIndexing)
	{
		for (unsigned int row=0; row<this->mColumns.size(); row++)
		{
			for (unsigned int col_count=0; col_count<this->mColumns[row].size(); col_count++)
			{
				rMatrix.AddValue(row, this->mColumns[row][col_count]-1,this->mValues[row][col_count]);
				if ((int)row!=this->mColumns[row][col_count]-1)
				{
				    rMatrix.AddValue(this->mColumns[row][col_count]-1, row, this->mValues[row][col_count]);
				}
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
                if ((int)row!=this->mColumns[row][col_count])
                {
                    rMatrix.AddValue(this->mColumns[row][col_count], row, this->mValues[row][col_count]);
                }
			}
		}
	}
}


//! @brief ... subtract two matrices
//! @param rOther ... Symmetric sparse matrix stored in the CSRVector2 format
//! @return reference to this matrix
template<class T>
NuTo::SparseMatrixCSRVector2Symmetric<T>& NuTo::SparseMatrixCSRVector2Symmetric<T>::operator-=  ( const NuTo::SparseMatrixCSRVector2Symmetric<T> &rOther )
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
			this->AddValue(row, rOther.mColumns[row][pos], - rOther.mValues[row][pos]);
		}
	}
	return *this;
}

//! @brief ... add two matrices
//! @param rOther ... Symmetric sparse matrix stored in the CSRVector2 format
//! @return reference to this matrix
template<class T>
NuTo::SparseMatrixCSRVector2Symmetric<T>& NuTo::SparseMatrixCSRVector2Symmetric<T>::operator+=  ( const NuTo::SparseMatrixCSRVector2Symmetric<T> &rOther )
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
			this->AddValue(row, rOther.mColumns[row][pos], rOther.mValues[row][pos]);
		}
	}
	return *this;
}

//! @brief ... matrix - matrix multiplication
//! @param rOther ... Symmetric sparse matrix stored in the CSR format
//! @return Symmetric sparse matrix stored in the CSR format
template<class T>
NuTo::SparseMatrixCSRVector2Symmetric<T> NuTo::SparseMatrixCSRVector2Symmetric<T>::operator* ( const NuTo::SparseMatrixCSRVector2Symmetric<T> &rOther ) const
{
    throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] To be implemented.");
/*  this is just copied from the general matrix (no symmtry)

    if (this->GetNumColumns() != rOther.GetNumRows())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");
	}
	SparseMatrixCSRVector2Symmetric<T> result(this->GetNumRows(), rOther.GetNumColumns());

	for (int thisRow = 0; thisRow < this->GetNumRows(); thisRow++)
	{
		const std::vector<T>& thisValueVec(this->mValues[thisRow]);
		const std::vector<int>& thisColumnVec(this->mColumns[thisRow]);
		for (unsigned int thisPos = 0; thisPos < thisValueVec.size(); thisPos++)
		{
			unsigned int thisColumn =thisColumnVec[thisPos];
			T thisValue = thisValueVec[thisPos];

			const std::vector<int>& otherColumnVec(rOther.mColumns[thisColumn]);
			const std::vector<T>& otherValueVec(rOther.mValues[thisColumn]);
			for (unsigned int otherPos = 0; otherPos < otherColumnVec.size(); otherPos++)
			{
				result.AddValue(thisRow, otherColumnVec[otherPos], thisValue * otherValueVec[otherPos]);
			}
		}
	}
	return result;
*/
}

//! @brief ... multiply sparse matrix with a full matrix
//! @param rFullMatrix ... full matrix which is multiplied with the sparse matrix
//! @return ... full matrix
template<class T>
NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> NuTo::SparseMatrixCSRVector2Symmetric<T>::operator* (const FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> &rMatrix) const
{
	if (this->GetNumColumns() != rMatrix.GetNumRows())
	{
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");
	}
	FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> result(this->GetNumRows(),rMatrix.GetNumColumns());
	result.setZero();
	if (this->HasOneBasedIndexing())
	{
		// loop over rows
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			// perform multiplication
			for (unsigned int pos = 0; pos < this->mColumns.size(); pos++)
			{
				int column = thisColumnVec[pos] - 1;
				T value = thisValueVec[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
				}
				if (column!=row)
				{
					//add the symmetric contribution
					for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
					{
						result(column,matrixCol) += value * rMatrix(row,matrixCol);
					}
				}
			}
		}
	}
	else
	{
		// loop over rows
		for (int row = 0; row < this->GetNumRows(); row++)
		{
			const std::vector<T>& thisValueVec(this->mValues[row]);
			const std::vector<int>& thisColumnVec(this->mColumns[row]);
			// perform multiplication
			for (unsigned int pos = 0; pos < thisColumnVec.size(); pos++)
			{
				int column = thisColumnVec[pos];
				T value = thisValueVec[pos];
				for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
				{
					result(row,matrixCol) += value * rMatrix(column,matrixCol);
					//std::cout << "add at " <<row << " " << matrixCol << " value " << value << " * " << rMatrix(column,matrixCol) << std::endl;
				}
				if (column!=row)
				{
					//add the symmetric contribution
					for (int matrixCol = 0; matrixCol < rMatrix.GetNumColumns(); matrixCol++)
					{
						result(column,matrixCol) += value * rMatrix(row,matrixCol);
						//std::cout << "add at symm " <<column << " " << matrixCol << " value " << value << " * " << rMatrix(row,matrixCol) << std::endl;
					}
				}
			}
		}
	}
	return result;
}

template<class T>
NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic> NuTo::SparseMatrixCSRVector2Symmetric<T>::TransMult(const NuTo::FullMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix) const
{
    throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] To be implemented.");
/*  this is just copied from the general matrix (no symmtry)
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
*/
}

//! @brief ... calculate the transpose of the matrix (transpose row and columns)
//! @return ... transpose of this matrix (sparse csr storage)
template<class T>
NuTo::SparseMatrixCSRVector2Symmetric<T> NuTo::SparseMatrixCSRVector2Symmetric<T>::Transpose() const
{
	return *this;
}

template<class T>
T NuTo::SparseMatrixCSRVector2Symmetric<T>::Sum() const
{
    T sum = 0;

    const auto& values = this->GetValues();
    const auto& columns = this->GetColumns();

    for (unsigned int row = 0; row < values.size(); row++)
    {
        for (unsigned int pos = 0; pos < values[row].size(); pos++)
        {
            sum += values[row][pos];

            unsigned int column = columns[row][pos];
            if (column != row)              // add off diagonals twice
                sum += values[row][pos];
        }
    }
    return sum;
}

//! @brief ... add the scaled other matrix
//! @param rOther ... other matrix
//! @param rFactor ... scalar factor
template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::AddScal(const SparseMatrixCSRVector2<T> &rOther, T rFactor)
{
    if (not rOther.IsSymmetric())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] rOther is not symmetric.");

    if ((this->GetNumColumns() != rOther.GetNumColumns()) || (this->GetNumRows() != rOther.GetNumRows()))
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid matrix dimensions.");

    if (this->HasOneBasedIndexing() || rOther.HasOneBasedIndexing())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] both matrices must have zero based indexing.");

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
void NuTo::SparseMatrixCSRVector2Symmetric<T>::AddScal(const SparseMatrixCSRSymmetric<T> &rOther, T rFactor)
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
    for (int row = 0; row < rOther.GetNumRows(); row++)
    {
        for (int pos = rOther.mRowIndex[row]; pos < rOther.mRowIndex[row + 1]; pos++)
        {
            this->AddValue(row, rOther.mColumns[pos], rFactor*rOther.mValues[pos]);
        }
    }

}

//! @brief ... reorder columns of the matrix
//! @param rMappingInitialToNewOrdering ... mapping fron initial to new ordering
template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering)
{
	throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] To be implemented.");
	/* copied from SparseMatrixCSR
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
*/
}

//! @brief ... resize matrix
//! @param rNumRows_ ... number of rows
template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::Resize(int rNumRows, int rNumColumns)
{
    // check for overflow
    assert(rNumRows < INT_MAX);
    assert(rNumRows >= 0);

    if (rNumRows!=rNumColumns)
		throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] number of rows and column has to be identical for symmetric matrices.");

    //no resize, since the reserved size is only decreased with a copy of a new object
    this->mValues  = std::vector<std::vector<T> >(rNumRows);
    this->mColumns = std::vector<std::vector<int > >(rNumRows);
}

//! @brief ... print info about the object
template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::Info() const
{
	std::cout << "number of rows: " << this->mValues.size() << std::endl;
	std::cout << "number of columns: " << this->mValues.size() << std::endl;
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
NuTo::SparseMatrixCSRVector2Symmetric<T> NuTo::SparseMatrixCSRVector2Symmetric<T>::Random(int rDimension, double rDensity, int rSeed)
{
    SparseMatrixCSRVector2General<T> matrix(rDimension, rDimension);

    // numValues(symmetric) is approx 0.5*numValues(general) since only one triangle is stored
    Matrix<T>::FillMatrixRandom(matrix, .5*rDensity, rSeed);
    return std::move(matrix.SymmetricPart());
}

template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::Add_TransA_B_C_Scal(
        const NuTo::SparseMatrixCSRVector2<T>& rA,
        const NuTo::SparseMatrixCSRVector2<T>& rB,
        const NuTo::SparseMatrixCSRVector2<T>& rC, T rScalar)
{

    // rA == rC, cheap asserts:
    assert(rA.GetNumColumns() == rC.GetNumColumns());
    assert(rA.GetNumRows()    == rC.GetNumRows());
    assert(rA.GetNumEntries() == rC.GetNumEntries());

    if (not rB.IsSymmetric())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] rB must be symmetric for the symmetric version of this operation");

    /*
     *  A.T * B * A
     *  where only the triangular matrix B§ is stored. Thus:
     *  B = B§ + B§.T - diag(B)
     *
     *     A.T * B  * A
     *  = [A.T * B$ * A] + [A.T * B$.T * A] - [A.T * diag(B) * A]
     *
     *  Because of the storage format, the following term is better
     *  = [B$.T * A].T * A + [B$ * A].T * A - [A.T * diag(B) * A]
     *
     *
     */




    int resultNumRows = rA.GetNumColumns();
    int resultNumCols = rA.GetNumColumns();

    if (rA.GetNumRows() != rB.GetNumRows())
    {
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong matrix dimensions for A.T * B");
    }
    if (resultNumRows != this->GetNumRows() or resultNumCols != this->GetNumColumns())
    {
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] matrix dimension mismatch of *this and A.T * B");
    }
    if (this->HasOneBasedIndexing() or rA.HasOneBasedIndexing() or rB.HasOneBasedIndexing())
    {
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] all matrices must have zero based indexing.");
    }



    // compressed:
    for (int bRow = 0; bRow < rB.GetNumRows(); ++bRow)
    {
        const auto& bValueVec = rB.GetValues()[bRow];
        const auto& bColumnVec = rB.GetColumns()[bRow];

        for (unsigned int iCol = 0; iCol < bColumnVec.size(); ++iCol)
        {
            int bCol = bColumnVec[iCol];
            T bValue = bValueVec[iCol];

            /*
             *  calculate
             *  [B$ * A].T * A
             */
            const auto& baColumnVec = rA.GetColumns()[bCol];
            const auto& baValueVec = rA.GetValues()[bCol];

            for (unsigned int i = 0; i < baColumnVec.size(); ++i)
            {
                int aCol = baColumnVec[i];
                // at this point: (B$ * A).T (aCol, bRow) = bValue * baValueVec[i] * rC


                const auto& aValueVec = rA.GetValues()[bRow];
                const auto& aColumnVec = rA.GetColumns()[bRow];

                for (unsigned int j = 0; j < aColumnVec.size(); ++j)
                    if (aCol <= aColumnVec[j]) // store symmetric part
                        this->AddValue(aCol, aColumnVec[j], bValue * baValueVec[i] * aValueVec[j] * rScalar);
            }

            /*
             *  skip the diagonals, equals
             *  - [A.T * diag(B) * A]
             */
            if (bCol == bRow)
                continue;

            /*
             *  calculate
             *  [B$.T * A].T * A
             */
            const auto& bTaColumnVec = rA.GetColumns()[bRow];
            const auto& bTaValueVec = rA.GetValues()[bRow];

            for (unsigned int i = 0; i < bTaColumnVec.size(); ++i)
            {
                int aCol = bTaColumnVec[i];
                // at this point: (B$.T * A).T (aCol, bCol) = bValue * bTaValueVec[i] * rC

                const auto& aValueVec = rA.GetValues()[bCol];
                const auto& aColumnVec = rA.GetColumns()[bCol];

                for (unsigned int j = 0; j < aColumnVec.size(); ++j)
                {
                    if (aCol <= aColumnVec[j]) // store symmetric part
                        this->AddValue(aCol, aColumnVec[j], bValue * bTaValueVec[i] * aValueVec[j] * rScalar);
                }
            }
        }
    }

}

template<class T>
void NuTo::SparseMatrixCSRVector2Symmetric<T>::Sub_TransA_B_Plus_C_D_Scal(
        const SparseMatrixCSRVector2<T>& rA,
        const SparseMatrixCSRVector2<T>& rB,
        const SparseMatrixCSRVector2<T>& rC,
        const SparseMatrixCSRVector2<T>& rD, T rScalar)
{
    if (
            rB.GetNumColumns() != rC.GetNumRows() or
            rC.GetNumColumns() != rB.GetNumRows() or
            rB.GetNumEntries() != rC.GetNumEntries())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] B.T should be equal to C for the symmetric version of this method.");

    if (
            rA.GetNumColumns() != rD.GetNumColumns() or
            rA.GetNumRows()    != rD.GetNumRows()    or
            rA.GetNumEntries() != rD.GetNumEntries())
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] A should be equal to D for the symmetric version of this method.");

    /*
     *  IGNORES B and D
     */

    /*
     *    A.T * C.T + C * A
     *  = (C * A).T + C * A
     *  can be generated with one B * A and clever AddValue
     */

    T minusScalar = -rScalar;

    int resultNumRows = rC.GetNumRows();
    int resultNumCols = rA.GetNumColumns();

    if (rC.GetNumColumns() != rA.GetNumRows())
    {
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] wrong matrix dimensions for A * B");
    }
    if (resultNumRows != this->GetNumRows() or resultNumCols != this->GetNumColumns())
    {
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] matrix dimension mismatch of *this and A * B");
    }
    if (this->HasOneBasedIndexing() or rA.HasOneBasedIndexing() or rB.HasOneBasedIndexing())
    {
        throw MathException(std::string("[") + __PRETTY_FUNCTION__ + "] all matrices must have zero based indexing.");
    }

    // calculate C * A
    for (int cRow = 0; cRow < rC.GetNumRows(); ++cRow)
    {
        const auto& cValueVec = rC.GetValues()[cRow];
        const auto& cColumnVec = rC.GetColumns()[cRow];

        for (unsigned int iCol = 0; iCol < cColumnVec.size(); ++iCol)
        {
            unsigned int bCol = cColumnVec[iCol];
            T bValue = cValueVec[iCol];

            const auto& aColumnVec = rA.GetColumns()[bCol];
            const auto& aValueVec = rA.GetValues()[bCol];

            for (unsigned int i = 0; i < aColumnVec.size(); ++i)
            {
                int aCol = aColumnVec[i];

                if (aCol == cRow)       // add diagonals twice
                    this->AddValue(aCol, cRow, bValue * aValueVec[i] * minusScalar * 2);

                else if (aCol < cRow)   // add off diagonals
                    this->AddValue(aCol, cRow, bValue * aValueVec[i] * minusScalar);

                else                    // add symmetric off diagonals
                    this->AddValue(cRow, aCol, bValue * aValueVec[i] * minusScalar);

            }
        }
    }

}


template <class T>
NuTo::SparseMatrixCSRVector2Symmetric<T>& NuTo::SparseMatrixCSRVector2Symmetric<T>::AsSparseMatrixCSRVector2Symmetric()
{
	return *this;
}

template <class T>
const NuTo::SparseMatrixCSRVector2Symmetric<T>& NuTo::SparseMatrixCSRVector2Symmetric<T>::AsSparseMatrixCSRVector2Symmetric() const
{
	return *this;
}

