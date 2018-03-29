#pragma once
#include "nuto/mechanics/dofs/DofInfo.h"
#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofMatrixSparse.h"

namespace NuTo
{

std::vector<DofType> AvailableDofTypes(const DofVector<int> n, std::vector<DofType> dofTypes)
{
    if (dofTypes.empty())
        return n.DofTypes();

    std::vector<DofType> intersection;
    for (const auto& dofTypeRow : n.DofTypes())
        for (const auto& dofTypeCol : dofTypes)
            if (dofTypeRow.Id() == dofTypeCol.Id())
            {
                intersection.push_back(dofTypeRow);
                continue;
            }
    return intersection;
}

class VectorAssembler
{
public:
    VectorAssembler(DofContainer<int> sizes = {})
    {
        Resize(sizes);
    }

    void Resize(DofContainer<int> sizes)
    {
        for (auto dofSize : sizes)
            mVector[dofSize.first].setZero(dofSize.second);
    }

    void Add(const DofVector<double>& v, const DofVector<int>& numbering, std::vector<DofType> dofTypes = {})
    {
        for (const auto& dofType : AvailableDofTypes(numbering, dofTypes))
        {
            for (unsigned localDofNumber = 0; localDofNumber < numbering[dofType].size(); ++localDofNumber)
            {
                const int globalDofNumber = numbering[dofType][localDofNumber];
                mVector[dofType][globalDofNumber] += v[dofType][localDofNumber];
            }
        }
    }

    void Reset()
    {
        for (const auto& dofType : mVector.DofTypes())
            mVector[dofType].setZero();
    }

    const DofVector<double>& Get() const
    {
        return mVector;
    }

private:
    DofVector<double> mVector;
};


class MatrixAssembler
{
public:
    MatrixAssembler(DofContainer<int> sizes = {})
    {
        Resize(sizes);
    }

    void Resize(DofContainer<int> sizes)
    {
        for (auto dofTypeRow : sizes)
            for (auto dofTypeCol : sizes)
                mMatrix(dofTypeRow.first, dofTypeCol.first).resize(dofTypeRow.second, dofTypeCol.second);
    }

    void Add(const DofMatrix<double>& m, const DofVector<int>& numbering, std::vector<DofType> dofTypes = {})
    {
        auto availableDofTypes = AvailableDofTypes(numbering, dofTypes);
        for (const auto& dofTypeRow : availableDofTypes)
        {
            for (const auto& dofTypeCol : availableDofTypes)
            {
                auto& tripletListDof = mTriplets(dofTypeRow, dofTypeCol);
                const auto& matrixDof = m(dofTypeRow, dofTypeCol);
                const auto& numberingRow = numbering[dofTypeRow];
                const auto& numberingCol = numbering[dofTypeCol];

                for (unsigned localRow = 0; localRow < numberingRow.size(); ++localRow)
                {
                    for (unsigned localCol = 0; localCol < numberingCol.size(); ++localCol)
                    {
                        const int globalRow = numberingRow[localRow];
                        const int globalCol = numberingCol[localCol];
                        const double value = matrixDof(localRow, localCol);
                        tripletListDof.push_back({globalRow, globalCol, value});
                    }
                }
            }
        }
    }

    static void Add(DofMatrixSparse<double>& rSparseMatrix, const DofMatrix<double>& m, const DofVector<int>& numbering,
                    std::vector<DofType> dofTypes = {})
    {
        auto availableDofTypes = AvailableDofTypes(numbering, dofTypes);
        for (const auto& dofTypeRow : availableDofTypes)
        {
            for (const auto& dofTypeCol : availableDofTypes)
            {
                auto& sparseMatrixDof = rSparseMatrix(dofTypeRow, dofTypeCol);
                const auto& matrixDof = m(dofTypeRow, dofTypeCol);
                const auto& numberingRow = numbering[dofTypeRow];
                const auto& numberingCol = numbering[dofTypeCol];

                for (unsigned localRow = 0; localRow < numberingRow.size(); ++localRow)
                {
                    for (unsigned localCol = 0; localCol < numberingCol.size(); ++localCol)
                    {
                        const int globalRow = numberingRow[localRow];
                        const int globalCol = numberingCol[localCol];
                        const double value = matrixDof(localRow, localCol);
                        sparseMatrixDof.coeffRef(globalRow, globalCol) += value;
                    }
                }
            }
        }
    }

    void Reset()
    {
        for (auto dofTypeRow : mMatrix.DofTypes())
            for (auto dofTypeCol : mMatrix.DofTypes())
                mMatrix(dofTypeRow, dofTypeCol).setZero();
    }

    void Finish()
    {
        for (auto dofTypeRow : mMatrix.DofTypes())
        {
            for (auto dofTypeCol : mMatrix.DofTypes())
            {
                const auto& tripletListDof = mTriplets(dofTypeRow, dofTypeCol);
                auto& matrix = mMatrix(dofTypeRow, dofTypeCol);
                matrix.setFromTriplets(tripletListDof.begin(), tripletListDof.end());
                matrix.makeCompressed();
            }
        }
    }

    const DofMatrixSparse<double>& Get() const
    {
        return mMatrix;
    }

private:
    using TripletList = std::vector<Eigen::Triplet<double>>;
    DofMatrixContainer<TripletList> mTriplets;

    DofMatrixSparse<double> mMatrix;
};


} /* NuTo */
