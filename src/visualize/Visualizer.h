#pragma once

#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "visualize/HandlerInterface.h"
#include "visualize/UnstructuredGrid.h"

#include "visualize/AverageHandler.h"
#include "visualize/VoronoiHandler.h"

namespace NuTo
{
namespace Visualize
{

//! Class to write visualization files of DOF values and cell data.
class Visualizer
{
public:
    //! Construct a visualizer with a group of cells to be visualized.
    //! @param cells Group of cells you want to visualize.
    //! @param handler implementation of the HandlerInterface
    Visualizer(Group<CellInterface>& cells, const HandlerInterface& handler)
        : mCells(cells)
        , mHandler(handler.Clone())
    {
        WriteGeometry();
    }

    //! Construct an average visualizer with a group of cells to be visualized.
    //! @param cells Group of cells you want to visualize.
    //! @param geometry average geometry
    Visualizer(Group<CellInterface>& cells, AverageGeometry geometry)
        : mCells(cells)
        , mHandler(std::make_unique<AverageHandler>(geometry))
    {
        WriteGeometry();
    }

    //! Construct an voronoi visualizer with a group of cells to be visualized.
    //! @param cells Group of cells you want to visualize.
    //! @param geometry voronoi geometry
    Visualizer(Group<CellInterface>& cells, VoronoiGeometry geometry)
        : mCells(cells)
        , mHandler(std::make_unique<VoronoiHandler>(geometry))
    {
        WriteGeometry();
    }


    //! Define DOF values that should be visualized.
    //! @param dof DofType to visualize.
    void DofValues(DofType dof)
    {
        int i = 0;
        for (const auto& cell : mCells)
        {
            mHandler->WriteDofValues(cell, dof, mPointIds[i], &mGrid);
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
            mHandler->CellData(i, values, name, &mGrid);
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
            mHandler->PointData(cell, f, mPointIds[i], name, &mGrid);
            ++i;
        }
    }

    //! Write out a VTK unstructured grid file.
    //! @param filename Name of the resulting file.
    //! @param asBinary ... true for output as binary vtu file
    void WriteVtuFile(std::string filename, bool asBinary = true)
    {
        mGrid.ExportVtuDataFile(filename, asBinary);
    }

private:
    Group<CellInterface>& mCells;
    std::vector<std::vector<int>> mPointIds;
    std::unique_ptr<HandlerInterface> mHandler;
    UnstructuredGrid mGrid;


    void WriteGeometry()
    {
        for (auto& cell : mCells)
            mPointIds.push_back(mHandler->WriteGeometry(cell, &mGrid));
    }
};

} // namespace Visualize
} // namespace NuTo
