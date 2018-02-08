#pragma once

#include <vector>
#include <memory>
#include <Eigen/Core>
#include "mechanics/dofs/DofType.h"

namespace NuTo
{

class CellInterface;

namespace Visualize
{

class UnstructuredGrid;

//! Interface for a handler class that visualizes a single cell.
class HandlerInterface
{
public:
    virtual std::unique_ptr<HandlerInterface> Clone() const = 0;

    //! Generate a visualize geometry for each cell and write it to the grid.
    //! @param cell Current cell to be visualized.
    //! @param grid Pointer to the grid where the geometry should be written to.
    virtual std::vector<int> WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid) = 0;

    //! Write DOF values into the grid.
    //! @param cell Current cell to be visualized.
    //! @param dof DofType to be visualized.
    //! @param pointIds IDs of the visualization points belonging to the cell.
    //! @param grid Pointer to the grid where the geometry should be written to.
    virtual void WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                UnstructuredGrid* grid) = 0;

    //! Write cell data into the grid.
    //! @param cellId ID of the cell to be visualized.
    //! @param values List of values to be visualized (e.g. integration point values).
    //! @param name Name to be used in the resulting output file for the data array.
    //! @param grid Pointer to the grid where the geometry should be written to.
    virtual void CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name,
                          UnstructuredGrid* grid) = 0;

    //! Write point data into the grid.
    //! @param cell Cell to be visualized.
    //! @param f Function to be visualized.
    //! @param pointIds IDs of the visualization points belonging to the cell.
    //! @param name Name to be used in the resulting output file for the data array.
    //! @param grid Pointer to the grid where the geometry should be written to.
    virtual void PointData(const CellInterface& cell, std::function<Eigen::VectorXd(Eigen::VectorXd)> f,
                           std::vector<int> pointIds, std::string name, UnstructuredGrid* grid) = 0;
};

} // namespace Visualize
} // namespace NuTo
