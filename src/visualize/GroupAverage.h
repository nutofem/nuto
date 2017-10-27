#pragma once

#include "base/Group.h"

#include "mechanics/cell/CellInterface.h"

#include "mechanics/elements/ElementCollection.h"

#include "visualize/CellGeometry.h"
#include "visualize/UnstructuredGrid.h"

namespace NuTo
{
namespace Visualize
{

using namespace NuTo::Groups;

class GroupAverage
{
public:
    GroupAverage(Group<CellInterface> cells, CellGeometry geometry)
        : mCells(cells)
        , mCellGeometry(geometry)
    {
    }

    // @brief defines the point and cell geometries of the averaged cells
    // @param We define loop through the cells and its geometry entries in the same fassion, as we do in Visualize
    // below. So there is no need to store ids in the UnstructuredGrid
    void ExtractGeometry()
    {
        for (auto& cell : mCells)
        {
            const auto& coordinateElement = cell.GetElementCollection().CoordinateElement();
            NodeValues coordinates = coordinateElement.ExtractNodeValues();

            std::vector<int> cornerIndices;
            for (auto naturalCornerCoordinate : mCellGeometry.mCornerCoords)
            {
                auto N = coordinateElement.GetNMatrix(naturalCornerCoordinate);
                int newPointId = mGrid.AddPoint(N * coordinates);
                cornerIndices.push_back(newPointId);
            }

            // the cell
            mGrid.AddCell(cornerIndices, mCellGeometry.mCellType);
        }
    }

    void Visualize(std::string file, Group<DofType> dofs, bool asBinary)
    {
        // register dof type names at the grid
        for (auto& dof : dofs)
            mGrid.DefineCellData(dof.GetName());

        for (auto& cell : mCells)
        {
        }
        mGrid.ExportVtuDataFile(file, asBinary);
    }

private:
    Group<CellInterface> mCells;
    CellGeometry mCellGeometry;
    UnstructuredGrid mGrid;
};
} /* Visualize */
} /* NuTo */
