// $Id$

#ifndef SPARSE_MATRIX_CSR_SYMMETRIC_H
#define SPARSE_MATRIX_CSR_SYMMETRIC_H

#include "nuto/math/Matrix.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
template <class T> class SparseMatrixCSRGeneral;
template <class T> class FullMatrix;
//! @author Stefan Eckardt, ISM
//! @date July 2009
//! @brief ... class for symmetric sparse matrices which are stored in CSR format
template <class T>
class SparseMatrixCSRSymmetric : public SparseMatrixCSR<T>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    //! @param rDimension_ ... dimension (number of rows and number of columns) of square matrix
    //! @param rNumReserveEntries_ ... number of entries for which memory is reserved (optional)
    SparseMatrixCSRSymmetric(int rDimension_ = 0, int rNumReserveEntries_= 0) : SparseMatrixCSR<T>(rDimension_, rNumReserveEntries_)
    {
    }

    //! @brief ... returns whether the matrix is symmetric or unsymmetric
    //! @return true if the matrix is symmetric and false if the matrix is unsymmetric
    bool IsSymmetric() const
    {
        return true;
    }

    //! @brief ... returns the number of columns
    //! @return number of columns
    int GetNumColumns() const
    {
        return this->mRowIndex.size() -  1;
    }

    //! @brief ... add nonzero entry to matrix
    //! @param rRow ... row of the nonzero entry (zero based indexing!!!)
    //! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
    //! @param rValue ... value of the nonzero entry
    void AddEntry(int rRow, int rColumn, T rValue)
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
    void Info() const
    {
        std::cout << "number of columns: " << this->mRowIndex.size() - 1 << std::endl;
        SparseMatrixCSR<T>::Info();
    }

    //! @brief ... import matrix from slang object stored in  a text file
    //! @param rFileName ... file name
    void ImportFromSLangText(const char* rFileName);

    //! @brief ... write nonzero matrix entries into a full matrix
    //! @param rFullMatrix ... the full matrix
    void WriteEntriesToFullMatrix(FullMatrix<T>& rFullMatrix) const
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

    //! @brief ... adds \f$(\boldsymbol{A}^T\,\boldsymbol{B}\,\boldsymbol{A})\f$ to the matrix
    //! @param rMatrixA ... matrix A (general sparse matrix in csr storage)
    //! @param rMatrixB ... matrix B (symmetric sparse matrix in csr storage)
    void Add_TransA_Mult_B_Mult_A(const NuTo::SparseMatrixCSRGeneral<T>& rMatrixA, const NuTo::SparseMatrixCSRSymmetric<T>& rMatrixB);

    //! @brief ... subtract t\f$(\boldsymbol{A}^T\boldsymbol{B}^T + \boldsymbol{B} \boldsymbol{A})\f$ from the matrix
    //! @param rMatrixA ... matrix A (general sparse matrix in csr storage)
    //! @param rMatrixB ... matrix B (general sparse matrix in csr storage)
    void Sub_TransA_Mult_TransB_Plus_B_Mult_A(const NuTo::SparseMatrixCSRGeneral<T>& rMatrixA, const NuTo::SparseMatrixCSRGeneral<T>& rMatrixB);

    //! @brief ... multiply sparse matrix with full matrix
    //! @param rMatrix ... full matrix
    //! @return ... result matrix (full storage)
    FullMatrix<T> operator* (const FullMatrix<T>& rMatrix) const;
    
    //! @brief ... multiplies the matrix with an scalar value
    //! @param rOther ... scalar value
    //! @return ... the multiplied matrix (sparse csr storage)
	SparseMatrixCSRSymmetric<T> operator* (const T& rScal) const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp("SparseMatrixCSR",boost::serialization::base_object< SparseMatrixCSR<T> >(*this));
    }

    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType)const;

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType);
#endif // ENABLE_SERIALIZATION
protected:
};
}
#endif // SPARSE_MATRIX_CSR_SYMMETRIC_H
