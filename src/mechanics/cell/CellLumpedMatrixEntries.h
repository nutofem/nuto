#pragma once

#include "mechanics/cell/DofVectorGenerator.h"
#include "mechanics/cell/CellInterface.h"
#include "base/Group.h"
#include <iostream>

namespace NuTo
{

class CellLumpedMatrixEntries : public DofVectorGenerator
{
public:
    CellLumpedMatrixEntries(Group<CellInterface> cells, CellInterface::MatrixFunction f, std::vector<DofType> dofs)
        : mCells(cells)
        , mF(f)
        , mDofs(dofs)
        , mCurrentCellIterator(mCells.begin())
    {
    }

    bool IsValid() const override
    {
        return mCurrentCellIterator != mCells.end();
    }

    void Next() override
    {
        mCurrentCellIterator++;
    }

    //! Transforms the local matrices calculated by mF into a diagonally lumped matrix
    //! @return vector entry that represents a local diagonal matrices
    //! @remark HRZ lumping is used here, the assumptions made are:
    //! - shape functions sum to 1, only then the total mass for a cell can be calculated by
    //!   summing over all entries of the local mass matrix
    //! - the total mass of a cell is the same for all components of the considered dof
    VectorEntry& Get() override
    {
        DofMatrix<double> localMatrix = mCurrentCellIterator->Integrate(mF);
        for (DofType dof : mDofs)
        {
            Eigen::VectorXd localDiagonalDof = localMatrix(dof, dof).diagonal();
            double diagonalSum = localDiagonalDof.sum();
            double fullSum = localMatrix(dof, dof).sum();
            localDiagonalDof *= (fullSum / diagonalSum) / dof.GetNum();

            mCurrentEntry.first[dof] = localDiagonalDof;
            mCurrentEntry.second[dof] = mCurrentCellIterator->DofNumbering(dof);
        }
        return mCurrentEntry;
    }

    std::vector<DofType> Dofs() const override
    {
        return mDofs;
    }

private:
    Group<CellInterface> mCells;
    CellInterface::MatrixFunction mF;
    std::vector<DofType> mDofs;

    Group<CellInterface>::ConstGroupIterator mCurrentCellIterator;
    VectorEntry mCurrentEntry;
};

} /* NuTo */
