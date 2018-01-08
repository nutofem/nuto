#pragma once

#include "mechanics/cell/CellInterface.h"
#include "visualize/UnstructuredGrid.h"

namespace NuTo
{
namespace Visualize
{

//! Cell handler that averages the cell values over each cell.
class QuadAverageHandler
{
public:
    //! Generate a visualize geometry for each cell and write it to the grid.
    //! @param cell Current cell to be visualized.
    //! @param grid Pointer to the grid where the geometry should be written to.
    std::vector<int> WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid);

    //! Write DOF values into the grid.
    //! @param cell Current cell to be visualized.
    //! @param dof DofType to be visualized.
    //! @param pointIds IDs of the visualization points belonging to the cell.
    //! @param grid Pointer to the grid where the geometry should be written to.
    void WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                        UnstructuredGrid* grid);

    //! Write cell data into the grid.
    //! @param cellId ID of the cell to be visualized.
    //! @param values List of values to be visualized (e.g. integration point values).
    //! @param name Name to be used in the resulting output file for the data array.
    //! @param grid Pointer to the grid where the geometry should be written to.
    void CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name, UnstructuredGrid* grid);

private:
    const std::vector<Eigen::VectorXd> cornerCoordinates = {
            Eigen::Vector3d(-1.0, -1.0, 0.0), Eigen::Vector3d(1.0, -1.0, 0.0), Eigen::Vector3d(1.0, 1.0, 0.0),
            Eigen::Vector3d(-1.0, 1.0, 0.0)};
};

} // namespace Visualize
} // namespace NuTo
