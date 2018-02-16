#pragma once

#include "mechanics/cell/DofVectorGenerator.h"
#include "mechanics/cell/CellInterface.h"
#include "base/Group.h"
#include <iostream>

namespace NuTo
{

class CellVectorEntries : public DofVectorGenerator
{
public:
    CellVectorEntries(const Group<CellInterface>& cells, CellInterface::VectorFunction f, std::vector<DofType> dofs)
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

    VectorEntry& Get() override
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
    const Group<CellInterface>& mCells;
    CellInterface::VectorFunction mF;
    std::vector<DofType> mDofs;

    Group<CellInterface>::ConstGroupIterator mCurrentCellIterator;
    VectorEntry mCurrentEntry;
};

} /* NuTo */
