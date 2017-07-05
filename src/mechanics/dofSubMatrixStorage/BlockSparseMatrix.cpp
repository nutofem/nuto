/*
 * BlockSparseMatrix.cpp
 *
 *  Created on: 6 Jan 2016
 *      Author: ttitsche
 */

#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "math/SparseMatrixCSRVector2General_Def.h"
#include "math/SparseMatrixCSRVector2Symmetric_Def.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/StructureBase.h"

#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrix.h"

#include "eigen3/Eigen/Sparse"


NuTo::BlockSparseMatrix::BlockSparseMatrix(const DofStatus& rDofStatus, bool rCanBeSymmetric) :
        BlockStorageBase(rDofStatus), mCanBeSymmetric(rCanBeSymmetric)
{
    AllocateSubmatrices();
}

//! @brief copy ctor
NuTo::BlockSparseMatrix::BlockSparseMatrix(const BlockSparseMatrix&  rOther) :
        BlockStorageBase(rOther.mDofStatus), mCanBeSymmetric(rOther.mCanBeSymmetric)
{
    for (const auto& it : rOther.mData)
    {
        if (it.second->IsSymmetric())
            mData[it.first] = std::unique_ptr<SparseMatrixCSRVector2<double>>(new SparseMatrixCSRVector2Symmetric<double>(it.second->AsSparseMatrixCSRVector2Symmetric()));
        else
            mData[it.first] = std::unique_ptr<SparseMatrixCSRVector2<double>>(new SparseMatrixCSRVector2General<double>(it.second->AsSparseMatrixCSRVector2General()));
    }
}

NuTo::BlockSparseMatrix::BlockSparseMatrix(NuTo::BlockSparseMatrix &&rOther)
    :BlockStorageBase(rOther.mDofStatus),
     mCanBeSymmetric(rOther.mCanBeSymmetric)
{
    mData =std::move(rOther.mData);
}

NuTo::BlockSparseMatrix::~BlockSparseMatrix()
{}

void NuTo::BlockSparseMatrix::AllocateSubmatrices()
{
    mData.clear();
    const auto& dofTypes = mDofStatus.GetDofTypes();
    for (auto dofRow : dofTypes)
        for (auto dofCol : dofTypes)
        {
            if (mCanBeSymmetric and dofRow == dofCol and mDofStatus.IsSymmetric(dofRow))
                mData[std::make_pair(dofRow, dofCol)] = std::unique_ptr<SparseMatrixCSRVector2<double>>(new SparseMatrixCSRVector2Symmetric<double>(0, 0));
            else
                mData[std::make_pair(dofRow, dofCol)] = std::unique_ptr<SparseMatrixCSRVector2<double>>(new SparseMatrixCSRVector2General<double>(0, 0));
            mData[std::make_pair(dofRow, dofCol)]->SetPositiveDefinite();
        }
}

void NuTo::BlockSparseMatrix::FixOffDiagonalDimensions()
{
    const auto& dofTypes = mDofStatus.GetDofTypes();
    for (auto dofRow : dofTypes)
        for (auto dofCol : dofTypes)
        {
            if (dofRow == dofCol)
                continue;
            if ((*this)(dofRow, dofCol).GetNumEntries() != 0)
                throw MechanicsException(__PRETTY_FUNCTION__, "You're about to resize a matrix with values inside. This should be wrong.");

            int numRows = (*this)(dofRow, dofRow).GetNumRows();
            int numCols = (*this)(dofCol, dofCol).GetNumColumns();
            (*this)(dofRow, dofCol).Resize(numRows, numCols);
        }
}

NuTo::BlockSparseMatrix& NuTo::BlockSparseMatrix::operator =(const BlockSparseMatrix& rOther)
{
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();
    for (auto dofRow : activeDofTypes)
        for (auto dofCol : activeDofTypes)
            (*this)(dofRow, dofCol) = rOther(dofRow, dofCol);
    return *this;
}

NuTo::BlockSparseMatrix& NuTo::BlockSparseMatrix::operator=(NuTo::BlockSparseMatrix &&rOther)
{
    mCanBeSymmetric = rOther.mCanBeSymmetric;
    mData = std::move(rOther.mData);
    return *this;
}

//NuTo::BlockSparseMatrix& NuTo::BlockSparseMatrix::operator =(BlockSparseMatrix&& rOther)
//{
//    const auto& activeDofs = mDofStatus.GetActiveDofTypes();
//    for (auto dofRow : activeDofs)
//        for (auto dofCol : activeDofs)
//            (*this)(dofRow, dofCol) = std::move(rOther(dofRow, dofCol));
//    return *this;
//}

NuTo::SparseMatrixCSRVector2<double>& NuTo::BlockSparseMatrix::operator ()(Node::eDof rDofRow, Node::eDof rDofCol)
{
    auto data = mData.find(std::make_pair(rDofRow, rDofCol));
    assert(data != mData.end());
    return *((*data).second);
}

const NuTo::SparseMatrixCSRVector2<double>& NuTo::BlockSparseMatrix::operator ()(Node::eDof rDofRow, Node::eDof rDofCol) const
{
    auto data = mData.find(std::make_pair(rDofRow, rDofCol));
    assert(data != mData.end());
    return *((*data).second);
}

NuTo::BlockFullVector<double> NuTo::BlockSparseMatrix::operator *(const BlockFullVector<double>& rRhs) const
{
    NuTo::BlockFullVector<double> result(mDofStatus);

    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();

    for (auto dofRow : activeDofTypes)
    {
        result[dofRow].resize((*this)(dofRow, dofRow).GetNumRows());
        result[dofRow].setZero();
        for (auto dofSum : activeDofTypes)
        {
            result[dofRow] += (*this)(dofRow, dofSum).operator*(rRhs[dofSum]);
            // Strange work-around: Calling the operator* directly instead of
            // using * fixes an error caused by eigen version 3.3~beta2-1
            // Eigen somehow tries to convert the SparseMatrixCSRVector2<double> to a double type and (obviously) fails.
        }
    }

    return result;
}

void NuTo::BlockSparseMatrix::AddScal(const BlockSparseMatrix& rRhs, double rScalar)
{
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();
    for (auto dofRow : activeDofTypes)
        for (auto dofCol : activeDofTypes)
        (*this)(dofRow, dofCol).AddScal(rRhs(dofRow, dofCol), rScalar);
}

void NuTo::BlockSparseMatrix::AddScalDiag(const BlockFullVector<double>& rRhs, double rScalar)
{
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();
    for (auto dof : activeDofTypes)
        (*this)(dof, dof).AddScalDiag(rRhs[dof], rScalar);
}

void NuTo::BlockSparseMatrix::CheckDimensions() const
{
    const auto& dofTypes = mDofStatus.GetDofTypes();
    /*
     * check for the same numRows in each row
     */
    for (auto dofRow : dofTypes)
    {
        auto numRowsReference = (*this)(dofRow, dofRow).GetNumRows();
        for (auto dofCol : dofTypes)
        {
            int numRows = (*this)(dofRow, dofCol).GetNumRows();
            if (numRows != numRowsReference)
            {
                std::stringstream s;
                s << "[" << __PRETTY_FUNCTION__ << "] Submatrix row dimension mismatch. \n";
                s << "(" << Node::DofToString(dofRow) << "," << Node::DofToString(dofCol) << ") has " << numRows << " rows \n";
                s << "(" << Node::DofToString(dofRow) << "," << Node::DofToString(dofRow) << ") has " << numRowsReference << " rows \n";
                throw MechanicsException(s.str());
            }
        }
    }

    /*
     * check for the same numCols in each column
     */
    for (auto dofCol : dofTypes)
    {
        auto numColsReference = (*this)(dofCol, dofCol).GetNumColumns();
        for (auto dofRow : dofTypes)
        {
            int numCols = (*this)(dofRow, dofCol).GetNumColumns();
            if (numCols != numColsReference)
            {
                std::stringstream s;
                s << "[" << __PRETTY_FUNCTION__ << "] Submatrix column dimension mismatch. \n";
                s << "(" << Node::DofToString(dofRow) << "," << Node::DofToString(dofCol) << ") has " << numCols << " columns \n";
                s << "(" << Node::DofToString(dofRow) << "," << Node::DofToString(dofRow) << ") has " << numColsReference << " columns \n";
                throw MechanicsException(s.str());
            }
        }
    }
}


int NuTo::BlockSparseMatrix::GetNumRowsDof(const std::set<Node::eDof>& rDofTypes) const
{
    int numRows = 0;
    for (auto dof : rDofTypes)
        numRows += (*this)(dof, dof).GetNumRows();

    return numRows;
}

int NuTo::BlockSparseMatrix::GetNumColumnsDof(const std::set<Node::eDof>& rDofTypes) const
{
    int numCols = 0;
    for (auto dof : rDofTypes)
        numCols += (*this)(dof, dof).GetNumColumns();

    return numCols;
}

void NuTo::BlockSparseMatrix::Info() const
{
    int minLength = 60;
    std::cout << "Num sub matrices: " << mData.size() << std::endl;
    for (auto& pair : mData)
    {
        std::string dofTypes = "[ " + Node::DofToString(pair.first.first) + " , " + Node::DofToString(pair.first.second) + " ]:";
        int numAdditionalBlanks = std::max(0, minLength - (int) dofTypes.length());
        const std::string& additionalBlanks = std::string(numAdditionalBlanks, ' ');
        const auto& matrix = pair.second;

        std::cout << dofTypes << additionalBlanks << "(" << (*matrix).GetNumRows() << "x" << (*matrix).GetNumColumns() << "): ";
        std::cout << "nonZeros: " << (*matrix).GetNumEntries() << std::endl;
    }
}

void NuTo::BlockSparseMatrix::SetZero()
{
    for (auto& pair : mData)
        pair.second->SetZeroEntries();
}

//! @brief inverts the matrix coefficient-wise
void NuTo::BlockSparseMatrix::CwiseInvert()
{
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();
    for (auto dofRow : activeDofTypes)
        for (auto dofCol : activeDofTypes)
            (*this)(dofRow, dofCol).CwiseInvert();
}


void NuTo::BlockSparseMatrix::Resize(const std::map<Node::eDof, int>& rNumRowDofsMap, const std::map<Node::eDof, int>& rNumColumnDofsMap)
{
    const auto& dofTypes = mDofStatus.GetDofTypes();
    for (auto dofRow : dofTypes)
        for (auto dofCol : dofTypes)
            (*this)(dofRow, dofCol).Resize(rNumRowDofsMap.at(dofRow), rNumColumnDofsMap.at(dofCol));
}

//! @brief adds the GetNumEntries() for all active dof types
int NuTo::BlockSparseMatrix::GetNumActiveEntires() const
{
    int numEntries = 0;
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();
    for (auto dofRow : activeDofTypes)
        for (auto dofCol : activeDofTypes)
        {
            auto data = mData.find(std::make_pair(dofRow, dofCol));
            if (data != mData.end())
            {
                numEntries += data->second->GetNumEntries();
            }
        }
    return numEntries;
}


bool NuTo::BlockSparseMatrix::HasSymmetricActiveDofs() const
{
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();

    for (auto dofRow : activeDofTypes)
        for (auto dofCol : activeDofTypes)
        {
            // diagonal blocks must be symmetric ...
            if (dofRow == dofCol and not (*this)(dofRow, dofCol).IsSymmetric())
                return false;

            // ... and off-diagonals must be zero
            if (dofRow != dofCol and (*this)(dofRow, dofCol).GetNumEntries() < 0)
                return false;
        }

    return true;
}

Eigen::MatrixXd NuTo::BlockSparseMatrix::ExportToFullMatrix() const
{
    CheckDimensions();
    const auto& activeDofTypes = mDofStatus.GetActiveDofTypes();
    Eigen::MatrixXd result (GetNumActiveRows(), GetNumActiveColumns());
    result.setZero();

    int blockStartRow = 0;
    for (auto dofRow : activeDofTypes)
    {
        int blockStartCol = 0;
        for (auto dofCol : activeDofTypes)
        {
            const auto& subMatrix = (*this)(dofRow, dofCol);
            const std::vector<std::vector<int> >& columns = subMatrix.GetColumns();
            const std::vector<std::vector<double>>& values = subMatrix.GetValues();

            // insert every nonzero element...
            for (unsigned int subMatrixRow = 0; subMatrixRow < columns.size(); ++subMatrixRow)
                for (unsigned int col_count=0; col_count < columns[subMatrixRow].size(); col_count++)
                {
                    int subMatrixCol = columns[subMatrixRow][col_count];
                    double subMatrixValue = values[subMatrixRow][col_count];

                    int exportRow = subMatrixRow + blockStartRow;
                    int exportCol = subMatrixCol + blockStartCol;

                    result(exportRow, exportCol) = subMatrixValue;
                    if (subMatrix.IsSymmetric() && exportRow != exportCol)
                        result(exportCol, exportRow) = subMatrixValue;
                }
            blockStartCol += subMatrix.GetNumColumns();
        }
        // every submatrix in the row dofRow has the same number of rows (CheckDimension)
        // --> one is picked and added to blockStartRow
        blockStartRow += (*this)(dofRow, dofRow).GetNumRows();
    }
    return result;
}

Eigen::SparseMatrix<double> NuTo::BlockSparseMatrix::ExportToEigenSparseMatrix() const
{

    SparseMatrixCSRGeneral<double> sparseMatrixCSRGeneral(ExportToCSRVector2General());

    std::vector<Eigen::Triplet<double>> tripletList;
    std::vector<double> val     = sparseMatrixCSRGeneral.GetValues();
    std::vector<int> colInd     = sparseMatrixCSRGeneral.GetColumns();
    std::vector<int> rowInd     = sparseMatrixCSRGeneral.GetRowIndex();

    for (unsigned i = 0; i < rowInd.size() - 1; ++i)
    {
        for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
            tripletList.push_back(Eigen::Triplet<double>(i, colInd[k], val[k]));
    }

    Eigen::SparseMatrix<double> eigenSparseMatrix(GetNumActiveRows(), GetNumActiveColumns());
    eigenSparseMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    eigenSparseMatrix.makeCompressed();

    return eigenSparseMatrix;
}


NuTo::SparseMatrixCSRVector2General<double> NuTo::BlockSparseMatrix::ExportToCSRVector2General() const
{
    SparseMatrixCSRVector2General<double> result;
    result.Resize(GetNumActiveRows(), GetNumActiveColumns());

    auto& vals = result.GetValuesReference();
    auto& cols = result.GetColumnsReference();

//    std::cout << GetNumActiveRows() << std::endl;
//    std::cout << vals.size()  << std::endl;

    int blockStartRow = 0;
    for (auto dofRow : mDofStatus.GetActiveDofTypes())
    {
        int numRowsDof = (*this)(dofRow, dofRow).GetNumRows();
        int blockStartCol = 0;
        // extract pointers to submatrices to avoid multiple access operations
        for (auto dofCol : mDofStatus.GetActiveDofTypes())
        {
            const auto& subMatrix = (*this)(dofRow, dofCol);
            const auto& dofVals = subMatrix.GetValues();
            const auto& dofCols = subMatrix.GetColumns();
//            std::cout << dofVals.size() << std::endl;

            if (subMatrix.IsSymmetric())
            {
                // insert non zero entries via AddValue
                for (unsigned int subMatrixRow = 0; subMatrixRow < dofCols.size(); ++subMatrixRow)
                    for (unsigned int col_count = 0; col_count < dofCols[subMatrixRow].size(); col_count++)
                    {
                        int subMatrixCol = dofCols[subMatrixRow][col_count];
                        double subMatrixValue = dofVals[subMatrixRow][col_count];

                        int exportRow = subMatrixRow + blockStartRow;
                        int exportCol = subMatrixCol + blockStartCol;

                        result.AddValue(exportRow, exportCol, subMatrixValue);      // normal value
                        if (exportRow != exportCol)                                 // exclude diagonals
                            result.AddValue(exportCol, exportRow, subMatrixValue);  // add transposed values
                    }
            }
            else
            {
                // concatenates the column/value vectors of each sub matrix
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                for (int iRow = 0; iRow < numRowsDof; ++iRow)
                {
                    int numEntries = vals[blockStartRow + iRow].size() + dofVals[iRow].size();
                    vals[blockStartRow + iRow].reserve(numEntries);
                    cols[blockStartRow + iRow].reserve(numEntries);
                    for (unsigned int iCol = 0; iCol < dofVals[iRow].size(); ++iCol)
                    {
                        vals[blockStartRow + iRow].push_back(dofVals[iRow][iCol]);
                        cols[blockStartRow + iRow].push_back(dofCols[iRow][iCol] + blockStartCol); // corrector for the column numbering
                    }
                }
            }
            blockStartCol += subMatrix.GetNumColumns();
        }
        blockStartRow += numRowsDof;
    }
    return result;
}

NuTo::SparseMatrixCSRGeneral<double> NuTo::BlockSparseMatrix::ExportToCSRGeneral() const
{
    return SparseMatrixCSRGeneral<double>(ExportToCSRVector2General());
}


std::unique_ptr<NuTo::SparseMatrixCSR<double>> NuTo::BlockSparseMatrix::ExportToCSR() const
{
    const auto& activeDofs = mDofStatus.GetActiveDofTypes();
    if (activeDofs.size() == 0)
        throw MechanicsException(__PRETTY_FUNCTION__, "No active dofs defined. Nothing to export.");

    auto dof = *activeDofs.begin();
    if (activeDofs.size() == 1 && mDofStatus.IsSymmetric(dof))
    {
        // symmetric case
        auto& ref = (*this)(dof,dof);
        assert(ref.IsSymmetric());
        std::unique_ptr<NuTo::SparseMatrixCSR<double>> matrix = std::make_unique<SparseMatrixCSRSymmetric<double>>(ref.AsSparseMatrixCSRVector2Symmetric());
        matrix->SetPositiveDefinite();
        return matrix;
    }
    else
    {
        return std::make_unique<SparseMatrixCSRGeneral<double>>(ExportToCSRVector2General());
    }
}


NuTo::SparseMatrixCSRVector2General<double> NuTo::BlockSparseMatrix::Get(std::string rDofRow, std::string rDofCol) const
{
    auto& ref = (*this)(Node::DofToEnum(rDofRow), Node::DofToEnum(rDofCol));
    if (ref.IsSymmetric())
        return ref.AsSparseMatrixCSRVector2Symmetric(); // calls appropriate Vector2General ctor
    else
        return ref.AsSparseMatrixCSRVector2General();
}


namespace NuTo
{
//! @brief output stream operator for outputs with cout or files
std::ostream& operator<<(std::ostream &rOut, const NuTo::BlockSparseMatrix &rBlockSparseMatrix)
{
    Eigen::IOFormat cleanFormat(Eigen::StreamPrecision, 0, " ", "\n", "|", " |");
    for (auto dof1 : rBlockSparseMatrix.mDofStatus.GetActiveDofTypes())
    {
        for (auto dof2 : rBlockSparseMatrix.mDofStatus.GetActiveDofTypes())
        {
            rOut << "[" << Node::DofToString(dof1) << " - " << Node::DofToString(dof2) << "]" << std::endl;
            rOut << rBlockSparseMatrix(dof1, dof2).ConvertToFullMatrix().format(cleanFormat) << std::endl;
        }
    }
    return rOut;
}
}
