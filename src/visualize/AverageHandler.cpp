#include "visualize/AverageHandler.h"

using namespace NuTo::Visualize;

AverageHandler::AverageHandler(AverageGeometry geometry)
    : mGeometry(geometry)
{
}

std::vector<int> AverageHandler::WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid)
{
    std::vector<int> pointIds;
    pointIds.reserve(mGeometry.cornerCoordinates.size());

    for (auto corner : mGeometry.cornerCoordinates)
        pointIds.push_back(grid->AddPoint(cell.Interpolate(corner)));
    grid->AddCell(pointIds, mGeometry.cellType);
    return pointIds;
}

void AverageHandler::WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                    UnstructuredGrid* grid)
{
    grid->DefinePointData(dof.GetName());
    for (size_t iCorner = 0; iCorner < mGeometry.cornerCoordinates.size(); ++iCorner)
        grid->SetPointData(pointIds[iCorner], dof.GetName(),
                           cell.Interpolate(mGeometry.cornerCoordinates[iCorner], dof));
}

void AverageHandler::CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name, UnstructuredGrid* grid)
{
    grid->DefineCellData(name);
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(values[0].rows(), values[0].cols());
    for (auto val : values)
        sum += val;

    grid->SetCellData(cellId, name, sum / values.size());
}
