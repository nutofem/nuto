#include <fstream>
#include <string>

#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "base/Exception.h"

namespace NuTo
{

// adds the product of trans(A) * B * A to the matrix (A is a general matrix, and B is a symmetric matrix)
template <>
void SparseMatrixCSRSymmetric<int>::Add_TransA_Mult_B_Mult_A(const NuTo::SparseMatrixCSRGeneral<int>&,
                                                             const NuTo::SparseMatrixCSRSymmetric<int>&)
{
    throw Exception("[SparseMatrixCSRSymmetric::Add_TransA_Mult_B_Mult_A] not implemented for this data-type.");
}

// adds the product of trans(A) * B * A to the matrix (A is a general matrix, and B is a symmetric matrix)
template <>
void SparseMatrixCSRSymmetric<double>::Add_TransA_Mult_B_Mult_A(const NuTo::SparseMatrixCSRGeneral<double>& rMatrixA,
                                                                const NuTo::SparseMatrixCSRSymmetric<double>& rMatrixB)
{
    assert(rMatrixB.GetNumRows() == rMatrixB.GetNumColumns());
    if (rMatrixB.GetNumRows() != rMatrixA.GetNumRows())
    {
        throw Exception(
                "[SparseMatrixCSRSymmetric::Add_TransA_Mult_B_Mult_A] invalid dimension of matrices A and B.");
    }
    assert(this->GetNumRows() == this->GetNumColumns());
    if (rMatrixA.GetNumColumns() != this->GetNumRows())
    {
        throw Exception(
                "[SparseMatrixCSRSymmetric::Add_TransA_Mult_B_Mult_A] invalid number of columns in matrix A.");
    }
    // loop over the columns of trans(A) == rows of A
    const std::vector<int>& matrixARowIndex = rMatrixA.GetRowIndex();
    const std::vector<int>& matrixAColumns = rMatrixA.GetColumns();
    const std::vector<double>& matrixAValues = rMatrixA.GetValues();
    for (int transAColumn = 0; transAColumn < rMatrixA.GetNumRows(); transAColumn++)
    {
        for (int transAPos = matrixARowIndex[transAColumn]; transAPos < matrixARowIndex[transAColumn + 1]; transAPos++)
        {
            int transARow = matrixAColumns[transAPos];
            double transAValue = matrixAValues[transAPos];

            // multiply with B matrix
            int BRow = transAColumn;

            // upper triangle and diagonal of B
            for (int BPos = rMatrixB.mRowIndex[BRow]; BPos < rMatrixB.mRowIndex[BRow + 1]; BPos++)
            {
                int BColumn = rMatrixB.mColumns[BPos];
                double BValue = rMatrixB.mValues[BPos];
                double transAValue_Mult_BValue = transAValue * BValue;
                // multiply A
                int ARow = BColumn;
                for (int APos = matrixARowIndex[ARow]; APos < matrixARowIndex[ARow + 1]; APos++)
                {
                    int AColumn = matrixAColumns[APos];
                    double AValue = matrixAValues[APos];
                    // add to this matrix (upper triangle storage)
                    int thisRow = transARow;
                    int thisColumn = AColumn;
                    if (thisColumn >= thisRow)
                    {
                        this->AddValue(thisRow, thisColumn, transAValue_Mult_BValue * AValue);
                    }
                }
            }

            // lower triangle of B
            for (int BColumn = 0; BColumn < transAColumn; BColumn++)
            {
                int BPos = rMatrixB.mRowIndex[BColumn];
                for (; BPos < rMatrixB.mRowIndex[BColumn + 1]; BPos++)
                {
                    if (rMatrixB.mColumns[BPos] >= transAColumn)
                    {
                        break;
                    }
                }
                if ((BPos < rMatrixB.mRowIndex[BColumn + 1]) && (rMatrixB.mColumns[BPos] == transAColumn))
                {
                    double BValue = rMatrixB.mValues[BPos];

                    double transAValue_Mult_BValue = transAValue * BValue;
                    // multiply A
                    int ARow = BColumn;
                    for (int APos = matrixARowIndex[ARow]; APos < matrixARowIndex[ARow + 1]; APos++)
                    {
                        int AColumn = matrixAColumns[APos];
                        double AValue = matrixAValues[APos];

                        // add to this matrix (upper triangle storage)
                        int thisRow = transARow;
                        int thisColumn = AColumn;
                        if (thisColumn >= thisRow)
                        {
                            this->AddValue(thisRow, thisColumn, transAValue_Mult_BValue * AValue);
                        }
                    }
                }
            }
        }
    }
}

// subtract (trans(A) * B + trans(B) * A) from the matrix (A and B are general matrices)
template <>
void SparseMatrixCSRSymmetric<int>::Sub_TransA_Mult_TransB_Plus_B_Mult_A(
        const NuTo::SparseMatrixCSRGeneral<int>&, const NuTo::SparseMatrixCSRGeneral<int>&)
{
    throw Exception(
            "[SparseMatrixCSRSymmetric::Sub_TransA_Mult_B_Plus_TransB_Mult_A] not implemented for this data-type.");
}

// subtract (trans(A) * trans(B) + B * A) from the matrix (A and B are general matrices)
template <>
void SparseMatrixCSRSymmetric<double>::Sub_TransA_Mult_TransB_Plus_B_Mult_A(
        const NuTo::SparseMatrixCSRGeneral<double>& rMatrixA, const NuTo::SparseMatrixCSRGeneral<double>& rMatrixB)
{
    if (rMatrixA.GetNumRows() != rMatrixB.GetNumColumns())
    {
        throw Exception("[SparseMatrixCSRSymmetric::Sub_TransA_Mult_TransB_Plus_B_Mult_A] invalid number of rows "
                            "of matrix A and B.");
    }
    if (rMatrixB.GetNumRows() != this->GetNumRows())
    {
        throw Exception("[SparseMatrixCSRSymmetric::Sub_TransA_Mult_TransB_Plus_B_Mult_A] invalid number of "
                            "columns of matrix B.");
    }
    if (rMatrixA.GetNumColumns() != this->GetNumColumns())
    {
        throw Exception("[SparseMatrixCSRSymmetric::Sub_TransA_Mult_TransB_Plus_B_Mult_A] invalid number of "
                            "columns of matrix A.");
    }

    // calculate B * A
    // matrix A
    const std::vector<int>& matrixARowIndex = rMatrixA.GetRowIndex();
    const std::vector<int>& matrixAColumns = rMatrixA.GetColumns();
    const std::vector<double> matrixAValues = rMatrixA.GetValues();

    // matrix B
    const std::vector<int>& matrixBRowIndex = rMatrixB.GetRowIndex();
    const std::vector<int>& matrixBColumns = rMatrixB.GetColumns();
    const std::vector<double> matrixBValues = rMatrixB.GetValues();

    for (int bRow = 0; bRow < rMatrixB.GetNumRows(); bRow++)
    {
        for (int bPos = matrixBRowIndex[bRow]; bPos < matrixBRowIndex[bRow + 1]; bPos++)
        {
            int bColumn = matrixBColumns[bPos];
            double bValue = matrixBValues[bPos];

            int aRow = bColumn;
            for (int aPos = matrixARowIndex[aRow]; aPos < matrixARowIndex[aRow + 1]; aPos++)
            {
                int aColumn = matrixAColumns[aPos];
                double aValue = matrixAValues[aPos];
                double value = -1.0 * aValue * bValue;
                if (aColumn >= bRow)
                {
                    // sub value from upper triangle
                    this->AddValue(bRow, aColumn, value);
                }
                if (bRow >= aColumn)
                {
                    // sub transpose from upper triangle
                    this->AddValue(aColumn, bRow, value);
                }
            }
        }
    }
}

// multiply sparse matrix with full matrix
template <>
Eigen::MatrixXi SparseMatrixCSRSymmetric<int>::operator*(const Eigen::MatrixXi&) const
{
    throw Exception("[SparseMatrixCSRSymmetric<int>::operator*] not implemented for this data type.");
}

// multiply sparse matrix with full matrix
template <>
Eigen::MatrixXd SparseMatrixCSRSymmetric<double>::operator*(const Eigen::MatrixXd& rMatrix) const
{
    if (this->GetNumColumns() != rMatrix.rows())
    {
        throw Exception("[SparseMatrixCSRSymmetric<int>::operator*] invalid number of rows in input matrix.");
    }
    Eigen::MatrixXd result(this->GetNumRows(), rMatrix.cols());

    for (int matrixCol = 0; matrixCol < rMatrix.cols(); matrixCol++)
    {
        const double* matrixValues = rMatrix.data() + matrixCol * rMatrix.rows();
        double* resultValues = result.data() + matrixCol * rMatrix.rows();
        for (int thisRow = 0; thisRow < this->GetNumRows(); thisRow++)
        {
            for (int thisPos = this->mRowIndex[thisRow] - mOneBasedIndexing;
                 thisPos < this->mRowIndex[thisRow + 1] - mOneBasedIndexing; thisPos++)
            {
                int thisColumn = this->mColumns[thisPos] - mOneBasedIndexing;
                double thisValue = this->mValues[thisPos];
                resultValues[thisRow] += thisValue * matrixValues[thisColumn];
                if (thisRow != thisColumn)
                {
                    resultValues[thisColumn] += thisValue * matrixValues[thisRow];
                }
            }
        }
    }
    return result;
}

// multiply sparse matrix with scalar
template <>
SparseMatrixCSRSymmetric<int> SparseMatrixCSRSymmetric<int>::operator*(const int&) const
{
    throw Exception("[SparseMatrixCSRSymmetric<int>::operator*] not implemented for this data type.");
}

// multiply sparse matrix with scalar
template <>
SparseMatrixCSRSymmetric<double> SparseMatrixCSRSymmetric<double>::operator*(const double& rScal) const
{
    SparseMatrixCSRSymmetric<double> result(*this);
    BOOST_FOREACH (double& val, result.mValues)
        val *= rScal;
    return result;
}

//! @brief ... add sparse matrix
//! @param rMatrix ... sparse matrix
//! @return ... this
template <>
SparseMatrix<int>& SparseMatrixCSRSymmetric<int>::operator+=(const SparseMatrixCSRSymmetric<int>&)
{
    throw Exception("[SparseMatrixCSRSymmetric<int>::operator+=] not implemented for this data type.");
}

//! @brief ... add sparse matrix
//! @param rMatrix ... sparse matrix
//! @return ... this
template <>
SparseMatrix<double>& SparseMatrixCSRSymmetric<double>::operator+=(const SparseMatrixCSRSymmetric<double>& rOther)
{
    if (this->GetNumColumns() != rOther.GetNumColumns())
    {
        throw Exception("[SparseMatrixCSRSymmetric<double>::operator*] invalid number of columns in input matrix.");
    }
    if (this->GetNumRows() != rOther.GetNumRows())
    {
        throw Exception("[SparseMatrixCSRSymmetric<double>::operator*] invalid number of rows in input matrix.");
    }


    for (int otherRow = 0; otherRow < rOther.GetNumRows(); otherRow++)
    {
        for (int otherPos = rOther.mRowIndex[otherRow] - mOneBasedIndexing;
             otherPos < rOther.mRowIndex[otherRow + 1] - mOneBasedIndexing; otherPos++)
        {
            this->AddValue(otherRow, rOther.mColumns[otherPos] - rOther.mOneBasedIndexing, rOther.mValues[otherPos]);
        }
    }

    return *this;
}

//! @brief ... add sparse matrix
//! @param rMatrix ... sparse matrix
//! @return ... this
template <>
SparseMatrix<int>& SparseMatrixCSRSymmetric<int>::operator+=(const SparseMatrixCSRVector2Symmetric<int>&)
{
    throw Exception("[SparseMatrixCSRSymmetric<int>::operator+=] not implemented for this data type.");
}


//! @brief ... add sparse matrix
//! @param rMatrix ... sparse matrix
//! @return ... this
template <>
SparseMatrix<double>& SparseMatrixCSRSymmetric<double>::
operator+=(const SparseMatrixCSRVector2Symmetric<double>& rOther)
{
    if (this->GetNumColumns() != rOther.GetNumColumns())
    {
        throw Exception("[SparseMatrixCSRSymmetric<double>::operator*] invalid number of columns in input matrix.");
    }
    if (this->GetNumRows() != rOther.GetNumRows())
    {
        throw Exception("[SparseMatrixCSRSymmetric<double>::operator*] invalid number of rows in input matrix.");
    }

    for (int otherRow = 0; otherRow < rOther.GetNumRows(); otherRow++)
    {
        for (unsigned int otherPos = 0; otherPos < rOther.mValues[otherRow].size(); otherPos++)
        {
            this->AddValue(otherRow, rOther.mColumns[otherRow][otherPos] - rOther.mOneBasedIndexing,
                           rOther.mValues[otherRow][otherPos]);
        }
    }

    return *this;
}


} // namespace NuTo
