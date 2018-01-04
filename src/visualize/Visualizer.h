#pragma once

#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "visualize/UnstructuredGrid.h"

namespace NuTo
{
namespace Visualize
{

template <typename THandler>
class Visualizer
{
public:
    Visualizer(Group<CellInterface>& cells)
        : mCells(cells)
    {
        for (auto& cell : mCells)
            mPointIds.push_back(mHandler.WriteGeometry(cell, &mGrid));
    }

    void DofValues(DofType dof)
    {
        int i = 0;
        for (const auto& cell : mCells)
        {
            mHandler.WriteDofValues(cell, dof, mPointIds[i], &mGrid);
            ++i;
        }
    }

    void CellData(std::function<Eigen::VectorXd(const CellData&, const CellIpData&)> f, std::string name)
    {
        int i = 0;
        for (const auto& cell : mCells)
        {
            auto values = cell.Eval(f);
            mHandler.CellData(i, values, name, &mGrid);
            ++i;
        }
    }

    void WriteVTKFile(std::string filename)
    {
        mGrid.ExportVtuDataFile(filename, false);
    }

private:
    Group<CellInterface>& mCells;
    std::vector<std::vector<int>> mPointIds;
    THandler mHandler;
    UnstructuredGrid mGrid;
};

} // namespace Visualize
} // namespace NuTo
