#pragma once

#include "base/ContainerView.h"
#include "mechanics/cell/CellInterface.h"
#include "visualize/UnstructuredGrid.h"

namespace NuTo
{
namespace Visualize
{

//! Class to write visualization files of DOF values and cell data.
//! @tparam THandler Cell handler class that converts computation cell into visualize output.
//!                  For examples, see QuadAverageHandler and TensorProductVoronoiHandler.
template <typename THandler>
class Visualizer
{
public:
    //! Construct a visualizer with a group of cells to be visualized.
    //! @param cells Group of cells you want to visualize.
    //! @param ...args Additional arguments that get passed on to the constructor of the cell handler.
    template <typename... HandlerArguments>
    Visualizer(ContainerView<CellInterface> cells, HandlerArguments&&... args)
        : mCells(cells)
        , mHandler(std::forward<HandlerArguments>(args)...)
    {
        for (auto& cell : mCells)
            mPointIds.push_back(mHandler.WriteGeometry(cell, &mGrid));
    }

    //! Define DOF values that should be visualized.
    //! @param dof DofType to visualize.
    void DofValues(DofType dof)
    {
        int i = 0;
        for (const auto& cell : mCells)
        {
            mHandler.WriteDofValues(cell, dof, mPointIds[i], &mGrid);
            ++i;
        }
    }

    //! Define cell data that should be visualized.
    //! @param f Function that is passed to the cell for evaluation.
    //! @param name Name to be used in the resulting output file for the data array.
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

    //! Visualize a function y = f(x) over a collection of cells
    //! @param f Function taking the coordinates as an Eigen vector and returning an Eigen vector
    //! @param name Name to be used in the resulting output file for the data array.
    void PointData(std::function<Eigen::VectorXd(Eigen::VectorXd)> f, std::string name)
    {
        int i = 0;
        for (const auto& cell : mCells)
        {
            mHandler.PointData(cell, f, mPointIds[i], name, &mGrid);
            ++i;
        }
    }

    //! Write out a VTK unstructured grid file.
    //! @param filename Name of the resulting file.
    void WriteVtuFile(std::string filename)
    {
        mGrid.ExportVtuDataFile(filename, false);
    }

private:
    ContainerView<CellInterface> mCells;
    std::vector<std::vector<int>> mPointIds;
    THandler mHandler;
    UnstructuredGrid mGrid;
};

} // namespace Visualize
} // namespace NuTo
