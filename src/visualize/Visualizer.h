#pragma once

#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "visualize/HandlerInterface.h"
#include "visualize/UnstructuredGrid.h"

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
    Visualizer(const Group<CellInterface>& cells, const HandlerInterface& handler);

    //! copy ctor that allows overriding an existing visualizer to e.g. change the cells
    Visualizer(Visualizer&&) = default;

    //! assigment operator that allows overriding an existing visualizer to e.g. change the cells
    Visualizer& operator=(Visualizer&&) = default;

    //! Define DOF values that should be visualized.
    //! @param dof DofType to visualize.
    void DofValues(DofType dof);

    //! Define cell data that should be visualized.
    //! @param f Function that is passed to the cell for evaluation.
    //! @param name Name to be used in the resulting output file for the data array.
    void CellData(std::function<Eigen::VectorXd(const CellIpData&)> f, std::string name);

    //! Visualize a function y = f(x) over a collection of cells
    //! @param f Function taking the coordinates as an Eigen vector and returning an Eigen vector
    //! @param name Name to be used in the resulting output file for the data array.
    void PointData(std::function<Eigen::VectorXd(Eigen::VectorXd)> f, std::string name);

    //! Write out a VTK unstructured grid file.
    //! @param filename Name of the resulting file.
    //! @param asBinary ... true for output as binary vtu file
    void WriteVtuFile(std::string filename, bool asBinary = true);

private:
    Group<CellInterface> mCells;
    std::vector<std::vector<int>> mPointIds;
    std::unique_ptr<HandlerInterface> mHandler;
    UnstructuredGrid mGrid;

    void WriteGeometry();
};
} // namespace Visualize
} // namespace NuTo
