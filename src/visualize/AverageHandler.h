#pragma once

#include "mechanics/cell/CellInterface.h"
#include "visualize/UnstructuredGrid.h"

namespace NuTo
{
namespace Visualize
{

//! Geometry description of one cell for the AverageHandler
struct AverageGeometry
{
    std::vector<Eigen::VectorXd> cornerCoordinates;
    eCellTypes cellType;
};

//! Cell handler that averages the cell values over each cell.
class AverageHandler
{
public:
    //! Construct for a given AverageGeometry
    //! @param geometry geometry description of one averaged cell
    AverageHandler(AverageGeometry geometry);

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

    //! Write point data into the grid.
    //! @param cell Cell to be visualized.
    //! @param f Function to be visualized.
    //! @param pointIds IDs of the visualization points belonging to the cell.
    //! @param name Name to be used in the resulting output file for the data array.
    //! @param grid Pointer to the grid where the geometry should be written to.
    void PointData(const CellInterface& cell, std::function<Eigen::VectorXd(Eigen::VectorXd)> f,
                   std::vector<int> pointIds, std::string name, UnstructuredGrid* grid);

private:
    AverageGeometry mGeometry;
};

} // namespace Visualize
} // namespace NuTo
