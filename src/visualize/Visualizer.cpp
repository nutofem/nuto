#include "visualize/Visualizer.h"

using namespace NuTo;
using namespace NuTo::Visualize;

Visualizer::Visualizer(const Group<CellInterface>& cells, const HandlerInterface& handler)
    : mCells(cells)
    , mHandler(handler.Clone())
{
    WriteGeometry();
}

void Visualizer::DofValues(DofType dof)
{
    int i = 0;
    for (const auto& cell : mCells)
    {
        mHandler->WriteDofValues(cell, dof, mPointIds[i], &mGrid);
        ++i;
    }
}

void Visualizer::CellData(std::function<Eigen::VectorXd(const class CellData&, const CellIpData&)> f, std::string name)
{
    int i = 0;
    for (const auto& cell : mCells)
    {
        auto values = cell.Eval(f);
        mHandler->CellData(i, values, name, &mGrid);
        ++i;
    }
}

void Visualizer::PointData(std::function<Eigen::VectorXd(Eigen::VectorXd)> f, std::string name)
{
    int i = 0;
    for (const auto& cell : mCells)
    {
        mHandler->PointData(cell, f, mPointIds[i], name, &mGrid);
        ++i;
    }
}

void Visualizer::WriteVtuFile(std::string filename, bool asBinary)
{
    mGrid.ExportVtuDataFile(filename, asBinary);
}

void Visualizer::WriteGeometry()
{
    for (auto& cell : mCells)
        mPointIds.push_back(mHandler->WriteGeometry(cell, &mGrid));
}
