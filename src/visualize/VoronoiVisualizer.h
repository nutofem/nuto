#pragma once
#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/elements/ElementCollection.h"
#include "visualize/CellGeometryVoronoi.h"
#include "visualize/UnstructuredGrid.h"

namespace NuTo
{
namespace Visualize
{

using namespace NuTo::Groups;

class GroupVoronoi
{
public:
    GroupVoronoi(Group<CellInterface> cells, CellGeometryVoronoi geometry)
        : mCells(cells)
        , mCellGeometry(geometry)
    {
    }

    void GeometryToGrid()
    {
        for (auto& cell : mCells)
        {
            const auto& coordinateElement = cell.GetElementCollection().CoordinateElement();

            // register all corner points of the cell
            for (auto incidenceCell : mCellGeometry.mVoronoiCells)
            {
                std::vector<int> cornerIndices;
                for (int id : incidenceCell)
                {
                    int newPointId =
                            mGrid.AddPoint(Interpolate(coordinateElement, mCellGeometry.mCellCornerCoords[id]));
                    cornerIndices.push_back(newPointId);
                }
                // register the cell
                mGrid.AddCell(cornerIndices, mCellGeometry.mCellType);
            }
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
            int cornerId;
            for (auto dof : dofs)
            {
                const auto& element = cell.GetElementCollection().DofElement(dof);
                cornerId = 0;
                for (auto incidenceCell : mCellGeometry.mVoronoiCells)
                {
                    for (int id : incidenceCell)
                    {
                        mGrid.SetPointData(currentPointId + cornerId, dof.GetName(),
                                           Interpolate(element, mCellGeometry.mCellCornerCoords[id]));
                        cornerId++;
                    }
                }
            }
            currentPointId += cornerId;
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

            for (unsigned int integrationPointID = 0; integrationPointID < vals.size();
                 ++integrationPointID, ++currentCellId)
            {
                for (unsigned int iType = 0; iType < vals[0].size(); ++iType)
                {
                    std::string name = vals[0][iType].name;
                    mGrid.SetCellData(currentCellId, name, vals[integrationPointID][iType].data);
                }
            }
        }
    }

    const UnstructuredGrid& GetUnstructuredGrid() const
    {
        return mGrid;
    }

private:
    Group<CellInterface> mCells;
    CellGeometryVoronoi mCellGeometry;
    UnstructuredGrid mGrid;
};

} /* Visualize */
} /* NuTo */
