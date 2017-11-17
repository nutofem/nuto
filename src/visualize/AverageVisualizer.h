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
    void GeometryToGrid()
    {
        for (auto& cell : mCells)
        {
            const auto& coordinateElement = cell.GetElementCollection().CoordinateElement();

            // register all corner points of the cell
            std::vector<int> cornerIndices;
            for (auto naturalCornerCoordinate : mCellGeometry.mCornerCoords)
            {
                int newPointId = mGrid.AddPoint(Interpolate(coordinateElement, naturalCornerCoordinate));
                cornerIndices.push_back(newPointId);
            }

            // register the cell
            mGrid.AddCell(cornerIndices, mCellGeometry.mCellType);
        }
    }

    void DofsToGrid(Group<DofType> dofs)
    {
        // register dof type names at the grid
        for (auto dof : dofs)
            mGrid.DefinePointData(dof.GetName());
        int currentPointId = 0;
        for (auto& cell : mCells)
        {
            for (auto dof : dofs)
            {
                const auto& element = cell.GetElementCollection().DofElement(dof);
                int cornerId = 0;
                for (auto naturalCornerCoordinate : mCellGeometry.mCornerCoords)
                {
                    mGrid.SetPointData(currentPointId + cornerId, dof.GetName(),
                                       Interpolate(element, naturalCornerCoordinate));
                    cornerId++;
                }
            }
            currentPointId += mCellGeometry.mCornerCoords.size();
        }
    }

    void IpValuesToGrid()
    {
        // register ip value names at the grid
        IpValues ipValuesForName = mCells.begin()->GetIpValues()[0];
        for (IpValue ipValue : ipValuesForName) // first cell, first integration pointipValues)
            mGrid.DefineCellData(ipValue.name);

        int currentCellId = 0;
        for (auto& cell : mCells)
        {
            // add average cell values
            std::vector<IpValues> vals = cell.GetIpValues();
            int numIps = vals.size();

            for (int iType = 0; iType < vals[0].size(); ++iType)
            {
                std::string name = vals[0][iType].name;
                Eigen::VectorXd data = vals[0][iType].data;
                for (int iIp = 1; iIp < numIps; ++iIp) // first one is in average
                    data += vals[iIp][iType].data;

                mGrid.SetCellData(currentCellId, name, data / numIps);
            }
            currentCellId++;
        }
    }

    const UnstructuredGrid& GetUnstructuredGrid() const
    {
        return mGrid;
    }

private:
    Group<CellInterface> mCells;
    CellGeometry mCellGeometry;
    UnstructuredGrid mGrid;
};
} /* Visualize */
} /* NuTo */
