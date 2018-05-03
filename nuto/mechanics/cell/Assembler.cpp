#include "Assembler.h"
#include "nuto/base/Exception.h"

std::vector<NuTo::DofType> AvailableDofTypes(const NuTo::DofVector<int>& n, std::vector<NuTo::DofType> dofTypes)
{
    if (dofTypes.empty())
        return n.DofTypes();

    std::vector<NuTo::DofType> intersection;
    for (const auto& dofTypeRow : n.DofTypes())
        for (const auto& dofTypeCol : dofTypes)
            if (dofTypeRow.Id() == dofTypeCol.Id())
            {
                intersection.push_back(dofTypeRow);
                continue;
            }
    return intersection;
}

void NuTo::Assembler::Add(DofVector<double>& rVector, const DofVector<double>& v, const DofVector<int>& numbering,
                          std::vector<DofType> dofTypes)
{
    for (const auto& dofType : AvailableDofTypes(numbering, dofTypes))
    {
        for (unsigned localDofNumber = 0; localDofNumber < numbering[dofType].size(); ++localDofNumber)
        {
            const int globalDofNumber = numbering[dofType][localDofNumber];
            rVector[dofType][globalDofNumber] += v[dofType][localDofNumber];
        }
    }
}

void NuTo::Assembler::Add(DofMatrixContainer<std::vector<Eigen::Triplet<double>>>& rTriplets,
                          const DofMatrix<double>& m, const DofVector<int>& numbering, std::vector<DofType> dofTypes)
{
    auto availableDofTypes = AvailableDofTypes(numbering, dofTypes);
    for (const auto& dofTypeRow : availableDofTypes)
    {
        for (const auto& dofTypeCol : availableDofTypes)
        {
            auto& tripletListDof = rTriplets(dofTypeRow, dofTypeCol);
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

void NuTo::Assembler::Add(DofMatrixSparse<double>& rSparseMatrix, const DofMatrix<double>& m,
                          const DofVector<int>& numbering, std::vector<DofType> dofTypes)
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


NuTo::VectorAssembler::VectorAssembler(DofContainer<int> sizes)
{
    Resize(sizes);
}

void NuTo::VectorAssembler::Resize(DofContainer<int> sizes)
{
    for (auto dofSize : sizes)
        mVector[dofSize.first].setZero(dofSize.second);
}

void NuTo::VectorAssembler::Add(const DofVector<double>& v, const DofVector<int>& numbering,
                                std::vector<DofType> dofTypes)
{
    Assembler::Add(mVector, v, numbering, dofTypes);
}

void NuTo::VectorAssembler::SetZero()
{
    for (const auto& dofType : mVector.DofTypes())
        mVector[dofType].setZero();
}

const NuTo::DofVector<double>& NuTo::VectorAssembler::Get() const
{
    return mVector;
}

NuTo::MatrixAssembler::MatrixAssembler(DofContainer<int> sizes)
{
    Resize(sizes);
}

void NuTo::MatrixAssembler::Resize(DofContainer<int> sizes)
{
    if (mFinished)
        throw Exception(__PRETTY_FUNCTION__,
                        "You are not allowed to resize the assembler matrix once ::Finished() is called.");

    for (auto dofTypeRow : sizes)
        for (auto dofTypeCol : sizes)
            mMatrix(dofTypeRow.first, dofTypeCol.first).resize(dofTypeRow.second, dofTypeCol.second);
}


void NuTo::MatrixAssembler::Add(const DofMatrix<double>& m, const DofVector<int>& numbering,
                                std::vector<DofType> dofTypes)
{
    if (mFinished)
        // We can now assemble in the exising nonzeros of this->mMatrix
        Assembler::Add(mMatrix, m, numbering, dofTypes);
    else
        Assembler::Add(mTriplets, m, numbering, dofTypes);
}


void NuTo::MatrixAssembler::SetZero()
{
    for (auto dofTypeRow : mMatrix.DofTypes())
        for (auto dofTypeCol : mMatrix.DofTypes())
            mMatrix(dofTypeRow, dofTypeCol).setZero();
}

void NuTo::MatrixAssembler::Finish()
{
    if (mFinished)
        return;

    for (auto dofTypeRow : mMatrix.DofTypes())
    {
        for (auto dofTypeCol : mMatrix.DofTypes())
        {
            auto& tripletListDof = mTriplets(dofTypeRow, dofTypeCol);
            auto& matrix = mMatrix(dofTypeRow, dofTypeCol);
            matrix.setFromTriplets(tripletListDof.begin(), tripletListDof.end());
            matrix.makeCompressed();
            tripletListDof.clear();
            tripletListDof.shrink_to_fit();
        }
    }
    mFinished = true;
}

const NuTo::DofMatrixSparse<double>& NuTo::MatrixAssembler::Get() const
{
    if (not mFinished)
        throw Exception(__PRETTY_FUNCTION__, "You have to call ::Finish() first to complete the assembly!");
    return mMatrix;
}
