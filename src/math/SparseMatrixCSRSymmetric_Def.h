#pragma once

#include "math/SparseMatrixCSR.h"

namespace NuTo
{
template <class T>
class SparseMatrixCSRGeneral;
template <class T>
class SparseMatrixCSRVector2Symmetric;
//! @author Stefan Eckardt, ISM
//! @date July 2009
//! @brief ... class for symmetric sparse matrices which are stored in CSR format
template <class T>
class SparseMatrixCSRSymmetric : public SparseMatrixCSR<T>
{
    friend class NuTo::SparseMatrixCSRVector2Symmetric<T>;

public:
    //! @brief ... constructor
    //! @param rDimension_ ... dimension (number of rows and number of columns) of square matrix
    //! @param rNumReserveEntries_ ... number of entries for which memory is reserved (optional)
    SparseMatrixCSRSymmetric(int rDimension = 0, int rNumReserveEntries = 0);

    //! @brief ... copy constructor
    SparseMatrixCSRSymmetric(const SparseMatrixCSRSymmetric<T>& rOther);

    //! @brief ... constructor
    SparseMatrixCSRSymmetric(const SparseMatrixCSRVector2Symmetric<T>& rOther);

    //! @brief ... returns whether the matrix is symmetric or unsymmetric
    //! @return true if the matrix is symmetric and false if the matrix is unsymmetric
    bool IsSymmetric() const override;

    //! @brief ... returns the number of columns
    //! @return number of columns
    int GetNumColumns() const override;

    void Info() const override;

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

    //! @brief ... adds \f$(\boldsymbol{A}^T\,\boldsymbol{B}\,\boldsymbol{A})\f$ to the matrix
    //! @param rMatrixA ... matrix A (general sparse matrix in csr storage)
    //! @param rMatrixB ... matrix B (symmetric sparse matrix in csr storage)
    void Add_TransA_Mult_B_Mult_A(const NuTo::SparseMatrixCSRGeneral<T>& rMatrixA,
                                  const NuTo::SparseMatrixCSRSymmetric<T>& rMatrixB);

    //! @brief ... subtract t\f$(\boldsymbol{A}^T\boldsymbol{B}^T + \boldsymbol{B} \boldsymbol{A})\f$ from the matrix
    //! @param rMatrixA ... matrix A (general sparse matrix in csr storage)
    //! @param rMatrixB ... matrix B (general sparse matrix in csr storage)
    void Sub_TransA_Mult_TransB_Plus_B_Mult_A(const NuTo::SparseMatrixCSRGeneral<T>& rMatrixA,
                                              const NuTo::SparseMatrixCSRGeneral<T>& rMatrixB);

    //! @brief ... multiply sparse matrix with full matrix
    //! @param rMatrix ... full matrix
    //! @return ... result matrix (full storage)
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    operator*(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix) const override;

    //! @brief ... multiplies the matrix with an scalar value
    //! @param rOther ... scalar value
    //! @return ... the multiplied matrix (sparse csr storage)
    SparseMatrixCSRSymmetric<T> operator*(const T& rScal) const;

    //! @brief ... add sparse matrix
    //! @param rMatrix ... sparse matrix
    //! @return ... this
    SparseMatrix<T>& operator+=(const SparseMatrixCSRVector2Symmetric<T>& rMatrix) override;

    //! @brief ... add sparse matrix
    //! @param rMatrix ... sparse matrix
    //! @return ... this
    SparseMatrix<T>& operator+=(const SparseMatrixCSRSymmetric<T>& rMatrix) override;

    //! @brief ... resize the matrix and initialize everything to zero
    //! @param  rRow ... number of rows
    //! @param  rCol ... number of columns
    void Resize(int rRow, int rCol) override;

    //! @brief ... resize the matrix and initialize everything to zero
    //! @param  rRow ... number of rows=number of columns
    void Resize(int rRow);

    SparseMatrixCSRSymmetric<T>& AsSparseMatrixCSRSymmetric() override;

#ifndef SWIG
    const SparseMatrixCSRSymmetric<T>& AsSparseMatrixCSRSymmetric() const override;
#endif

protected:
};
}
