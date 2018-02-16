#pragma once

#include "mechanics/cell/DofMatrixGenerator.h"
#include "mechanics/cell/CellInterface.h"
#include "base/Group.h"
#include <iostream>

namespace NuTo
{

class CellMatrixEntries : public DofMatrixGenerator
{
public:
    CellMatrixEntries(Group<CellInterface> cells, CellInterface::MatrixFunction f, std::vector<DofType> dofs)
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

    MatrixEntry& Get() override
    {
        mCurrentEntry.first = mCurrentCellIterator->Integrate(mF);
        for (DofType dof : mDofs)
            mCurrentEntry.second[dof] = mCurrentCellIterator->DofNumbering(dof);
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
    MatrixEntry mCurrentEntry;
};

} /* NuTo */
