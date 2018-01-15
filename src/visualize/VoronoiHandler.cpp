#include "visualize/VoronoiHandler.h"

using namespace NuTo::Visualize;

VoronoiHandler::VoronoiHandler(VoronoiGeometry geometry)
    : mGeometry(geometry)
{
}

std::vector<int> VoronoiHandler::WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid)
{
    std::vector<int> pointIds;
    for (auto pointCoords : mGeometry.pointCoordinates)
        pointIds.push_back(grid->AddPoint(cell.Interpolate(pointCoords)));

    std::vector<int> localSubCells;
    for (auto voronoiCell : mGeometry.voronoiCells)
    {
        std::vector<int> localIds;
        for (int i : voronoiCell.cellCornerIds)
            localIds.push_back(pointIds[i]);
        localSubCells.push_back(grid->AddCell(localIds, voronoiCell.cellType));
    }
    mSubCells.push_back(localSubCells);
    return pointIds;
}

void VoronoiHandler::WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                    UnstructuredGrid* grid)
{
    grid->DefinePointData(dof.GetName());
    for (size_t iPoint = 0; iPoint < mGeometry.pointCoordinates.size(); ++iPoint)
        grid->SetPointData(pointIds[iPoint], dof.GetName(), cell.Interpolate(mGeometry.pointCoordinates[iPoint], dof));
}

void VoronoiHandler::CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name, UnstructuredGrid* grid)
{
    assert(values.size() == mSubCells[cellId].size());
    grid->DefineCellData(name);
    for (size_t iValue = 0; iValue < values.size(); ++iValue)
        grid->SetCellData(mSubCells[cellId][iValue], name, values[iValue]);
}
