// $Id: SparseMatrixCSRVector2General_Def.h 235 2010-04-22 09:25:38Z arnold2 $

#ifndef SPARSE_MATRIX_CSR_VECTOR2_GENERAL_DEF_H
#define SPARSE_MATRIX_CSR_VECTOR2_GENERAL_DEF_H

#include "nuto/math/SparseMatrixCSRVector2.h"

//the template class SparseMatrixCSRVector2General is split into definition and implementation, include only the definition here
#include "nuto/math/FullMatrix.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
template <class T> class SparseMatrixCSRGeneral;
//! @author Stefan Eckardt, ISM
//! @date July 2009
//! @brief ... class for general sparse matrices which are stored in CSR format
//!            this file contains the function definitions, whereas the implementation goes either into SparseMatrixCSRVector2General.h or the cpp file
//!            this was necessary to avoid cyclic inclusions of header templates

template <class T>
class SparseMatrixCSRVector2General : public SparseMatrixCSRVector2<T>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    //! @param rNumRows_ ... number of rows
    //! @param rNumColumns_ ... number of columns
    SparseMatrixCSRVector2General(int rNumRows_=0, int rNumColumns_=0);

    //! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a predefined tolerance)
    //! @param rFullMatrix ... input matrix (full storage)
    //! @param rAbsoluteTolerance ... absolute tolerance
    //! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance * max(abs(rMatrixEntry))
    SparseMatrixCSRVector2General(FullMatrix<T>& rFullMatrix, double rAbsoluteTolerance = 0, double rRelativeTolerance = 1e-14);

    //! @brief ... create sparse matrix with vector of vector from standard CSR format
    //! @param rCSRMatrix ... input matrix (full storage)
    SparseMatrixCSRVector2General(const SparseMatrixCSRGeneral<T>& rCSRMatrix);

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const;

    //! @brief ... returns whether the matrix is symmetric or unsymmetric
    //! @return true if the matrix is symmetric and false if the matrix is unsymmetric
    bool IsSymmetric() const;

    //! @brief ... add nonzero entry to matrix
    //! @param rRow ... row of the nonzero entry (zero based indexing!!!)
    //! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
    //! @param rValue ... value of the nonzero entry
    void AddEntry(int rRow, int rColumn, T rValue);

    //! @brief ... import matrix from slang object stored in  a text file
    //! @param rFileName ... file name
    void ImportFromSLangText(const char* rFileName);

    //! @brief ... write nonzero matrix entries into a full matrix
    //! @param rFullMatrix ... the full matrix
    void WriteEntriesToFullMatrix(FullMatrix<T>& rFullMatrix) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param rType ... type of file, either BINARY, XML or TEXT
    void Save ( const std::string &filename, std::string rType)const;

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore ( const std::string &filename,  std::string rType);
#endif // ENABLE_SERIALIZATION


    //! @brief ... add two matrices
    //! @param rOther ... general sparse matrix stored in the CSRVector2 format
    //! @return general sparse matrix stored in the CSR format
    SparseMatrixCSRVector2General<T> operator+ ( const SparseMatrixCSRVector2General<T> &rOther );

    //! @brief ... subtract two matrices
    //! @param rOther ... general sparse matrix stored in the CSRVector2 format
    //! @return general sparse matrix stored in the CSR format
    SparseMatrixCSRVector2General<T> operator- ( const SparseMatrixCSRVector2General<T> &rOther );

    //! @brief ... subtract two matrices
    //! @param rOther ... general sparse matrix stored in the CSRVector2 format
    //! @return reference to this matrix
    SparseMatrixCSRVector2General<T>& operator-=  ( const SparseMatrixCSRVector2General<T> &rOther );

    //! @brief ... add two matrices
    //! @param rOther ... general sparse matrix stored in the CSRVector2 format
    //! @return reference to this matrix
    SparseMatrixCSRVector2General<T>& operator+=  ( const SparseMatrixCSRVector2General<T> &rOther );

    //! @brief ... matrix - matrix multiplication
    //! @param rOther ... general sparse matrix stored in the CSR format
    //! @return general sparse matrix stored in the CSR format
    SparseMatrixCSRVector2General<T> operator* ( const SparseMatrixCSRVector2General<T> &rOther ) const;

    //! @brief ... multiplies the matrix with an scalar value
    //! @param rOther ... scalar value
    //! @return ... the multiplied matrix (sparse csr storage)
    SparseMatrixCSRVector2General<T> operator* ( const T &rOther ) const;

    //! @brief ... multiply sparse matrix with a full matrix
    //! @param rFullMatrix ... full matrix which is multiplied with the sparse matrix
    //! @return ... full matrix
    FullMatrix<T> operator* (const FullMatrix<T> &rMatrix) const;

    FullMatrix<T> TransMult(const FullMatrix<T>& rMatrix) const;

    //! @brief ... calculate the transpose of the matrix (transpose row and columns)
    //! @return ... transpose of this matrix (sparse csr storage)
    SparseMatrixCSRVector2General<T> Transpose() const;

    //! @brief ... perform Gauss algorithm (matrix and right hand side are reordered and modified)
    //! @param rRhs ... right-hand side vector (input and output object)
    //! @param rMappingNewToInitialOrdering ... mapping from new ordering to initial ordering (output object)
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to new ordering (output object)
    //! @param rRelativeTolerance ... relative tolerance for zero matrix entries
    void Gauss(FullMatrix<T>& rRhs, std::vector<int>& rMappingNewToInitialOrdering, std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance = 1e-14);

    //! @brief ... reorder columns of the matrix
    //! @param rMappingInitialToNewOrdering ... mapping fron initial to new ordering
    void ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering);

    //! @brief ... Concatenate columns from another matrix
    //! @param rOther ... other matrix with same number of rows
    void ConcatenateColumns(const SparseMatrixCSRVector2General<T>& rOther);

    //! @brief ... Concatenate rows from another matrix
    //! @param rOther ... other matrix with same number of columns
    void ConcatenateRows(const SparseMatrixCSRVector2General<T>& rOther);
protected:
};
}
#endif // SPARSE_MATRIX_CSR_VECTOR2_GENERAL_DEF_H
