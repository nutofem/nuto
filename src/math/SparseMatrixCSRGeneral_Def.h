#pragma once

#include "math/SparseMatrixCSR.h"

namespace NuTo
{
template <class T>
class SparseMatrixCSRVector2General;

//! @author Stefan Eckardt, ISM
//! @date July 2009
//! @brief ... class for general sparse matrices which are stored in CSR format
template <class T>
class SparseMatrixCSRGeneral : public SparseMatrixCSR<T>
{
    friend class NuTo::SparseMatrixCSRVector2General<T>;

public:
    //! @brief ... constructor
    //! @param rNumRows_ ... number of rows
    //! @param rNumColumns_ ... number of columns
    //! @param rNumReserveEntries_ ... number of entries for which memory is reserved (optional)
    SparseMatrixCSRGeneral(int rNumRows_ = 0, int rNumColumns_ = 0, unsigned int rNumReserveEntries_ = 0);

    //! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a
    //! predefined tolerance)
    //! @param rFullMatrix ... input matrix (full storage)
    //! @param rAbsoluteTolerance ... absolute tolerance
    //! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance *
    //! max(abs(rMatrixEntry))
    SparseMatrixCSRGeneral(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rFullMatrix,
                           double rAbsoluteTolerance = 0, double rRelativeTolerance = 1e-14);

    //! @brief ... create sparse matrix from full matrix (considers only matrix entries which absolute value exceeds a
    //! predefined tolerance)
    //! @param rFullMatrix ... input matrix (full storage)
    //! @param rAbsoluteTolerance ... absolute tolerance
    //! @param rRelative tolerance ... relative tolerance (tolerance = rAbsoluteTolerance + rRelativeTolerance *
    //! max(abs(rMatrixEntry))
    SparseMatrixCSRGeneral(const NuTo::SparseMatrixCSRVector2General<T>& rCSR2Matrix);

    void Info() const override;

    //! @brief ... resize matrix
    //! @param rNumRows_ ... number of rows
    //! @param rNumColumns_ ... number of columns
    void Resize(int rNumRows_, int rNumColumns_) override;

    //! @brief ... returns whether the matrix is symmetric or unsymmetric
    //! @return true if the matrix is symmetric and false if the matrix is unsymmetric
    bool IsSymmetric() const override;

    //! @brief ... returns the number of columns
    //! @return number of columns
    int GetNumColumns() const override;

    //! @brief ... add nonzero entry to matrix
    //! @param rRow ... row of the nonzero entry (zero based indexing!!!)
    //! @param rColumn ... column of the nonzero entry (zero based indexing!!!)
    //! @param rValue ... value of the nonzero entry
    void AddValue(int rRow, int rColumn, const T& rValue) override;

    //! @brief ... return the matrix type
    NuTo::eSparseMatrixType GetSparseMatrixType() const override;

    //! @brief ... write nonzero matrix entries into a matrix
    //! @return ... the matrix
    virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ConvertToFullMatrix() const override;

    //! @brief ... add two matrices
    //! @param rOther ... general sparse matrix stored in the CSR format
    //! @return general sparse matrix stored in the CSR format
    SparseMatrixCSRGeneral<T> operator+(const NuTo::SparseMatrixCSRGeneral<T>& rOther) const;

    //! @brief ... subtract two matrices
    //! @param rOther ... general sparse matrix stored in the CSR format
    //! @return general sparse matrix stored in the CSR format
    SparseMatrixCSRGeneral<T> operator-(const NuTo::SparseMatrixCSRGeneral<T>& rOther) const;

    //! @brief ... subtract two matrices
    //! @param rOther ... general sparse matrix stored in the CSR format
    //! @return reference to this matrix
    SparseMatrixCSRGeneral<T>& operator-=(const NuTo::SparseMatrixCSRGeneral<T>& rOther);


    //! @brief ... add two matrices
    //! @param rOther ... general sparse matrix stored in the CSR format
    //! @return reference to this matrix
    SparseMatrixCSRGeneral<T>& operator+=(const NuTo::SparseMatrixCSRGeneral<T>& rOther);

    //! @brief ... matrix - matrix multiplication
    //! @param rOther ... general sparse matrix stored in the CSR format
    //! @return general sparse matrix stored in the CSR format
    SparseMatrixCSRGeneral<T> operator*(const NuTo::SparseMatrixCSRGeneral<T>& rOther) const;

    //! @brief ... multiplies the matrix with an scalar value
    //! @param rOther ... scalar value
    //! @return ... the multiplied matrix (sparse csr storage)
    SparseMatrixCSRGeneral<T> operator*(const T& rOther) const;

    //! @brief ... multiply sparse matrix with a full matrix
    //! @param rFullMatrix ... full matrix which is multiplied with the sparse matrix
    //! @return ... full matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    operator*(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix) const override;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    TransMult(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix) const;

    //! @brief ... calculate the transpose of the matrix (transpose row and columns)
    //! @return ... transpose of this matrix (sparse csr storage)
    SparseMatrixCSRGeneral<T> Transpose() const;

    //! @brief ... remove entry from matrix
    //! @param rRow ... row index (zero based indexing!!!)
    //! @param rColumn ... column index (zero based indexing!!!)
    void RemoveEntry(int rRow, int rColumn);

    //! @brief ... remove columns from the end of the matrix
    //! @param rNumColumn ... number of colums to be removed
    void RemoveLastColumns(unsigned int rNumColumns);

    //! @brief ... perform Gauss algorithm (matrix and right hand side are reordered and modified)
    //! @param rRhs ... right-hand side vector (input and output object)
    //! @param rMappingNewToInitialOrdering ... mapping from new ordering to initial ordering (output object)
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to new ordering (output object)
    //! @param rRelativeTolerance ... relative tolerance for zero matrix entries
    void Gauss(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rRhs, std::vector<int>& rMappingNewToInitialOrdering,
               std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance = 1e-14);

    //! @brief ... perform Gauss algorithm (matrix and right hand side are reordered and modified)
    //! @param rRhs ... right-hand side matrix (input and output object, kind of multiple rhs)
    //! @param rMappingNewToInitialOrdering ... mapping from new ordering to initial ordering (output object)
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to new ordering (output object)
    //! @param rRelativeTolerance ... relative tolerance for zero matrix entries
    void Gauss(NuTo::SparseMatrixCSRGeneral<T>& rRhs, std::vector<int>& rMappingNewToInitialOrdering,
               std::vector<int>& rMappingInitialToNewOrdering, double rRelativeTolerance = 1e-14);

    //! @brief ... compute the maximum eigenvalue and eigenvector
    //! @param rStart ... starting vector for the iteration should be ||rStart|| = 1
    //! @param tol ... relative error for the iteration
    //! @return the maximum eigenvalue, rStart will become the eigenvector to the maximum eigenvalue
    void GetMaximumEigenvalueAndEigenvector(Eigen::Matrix<T, Eigen::Dynamic, 1>& rStart, T& maximumEigenvalue,
                                            double tol = 1.e-6);

    //! @brief ... reorder columns of the matrix
    //! @param rMappingInitialToNewOrdering ... mapping fron initial to new ordering
    void ReorderColumns(const std::vector<int>& rMappingInitialToNewOrdering);

    SparseMatrixCSRGeneral<T>& AsSparseMatrixCSRGeneral() override;
#ifndef SWIG
    const SparseMatrixCSRGeneral<T>& AsSparseMatrixCSRGeneral() const override;
#endif


protected:
    //! @brief ... number of columns
    int mNumColumns = 0;
};
}
