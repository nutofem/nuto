#include "QuadAverageHandler.h"

using namespace NuTo::Visualize;

std::vector<int> QuadAverageHandler::WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid)
{
    std::vector<int> pointIds;
    for (auto corner : cornerCoordinates)
        pointIds.push_back(grid->AddPoint(cell.Interpolate(corner)));
    grid->AddCell(pointIds, NuTo::eCellTypes::QUAD);
    return pointIds;
}

void QuadAverageHandler::WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                        UnstructuredGrid* grid)
{
    grid->DefinePointData(dof.GetName());
    int i = 0;
    for (auto corner : cornerCoordinates)
    {
        grid->SetPointData(pointIds[i], dof.GetName(), cell.Interpolate(corner, dof));
        ++i;
    }
}

void QuadAverageHandler::CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name,
                                  UnstructuredGrid* grid)
{
    grid->DefineCellData(name);
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(values[0].rows(), values[0].cols());
    for (auto val : values)
    {
        sum += val;
    }
    grid->SetCellData(cellId, name, sum / values.size());
}
