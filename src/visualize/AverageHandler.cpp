#include "visualize/AverageHandler.h"
#include "visualize/UnstructuredGrid.h"
#include "mechanics/cell/CellInterface.h"
#include "math/EigenCompanion.h"

using namespace NuTo::Visualize;

AverageHandler::AverageHandler(AverageGeometry geometry)
    : mGeometry(geometry)
{
}

std::unique_ptr<HandlerInterface> AverageHandler::Clone() const
{
    return std::make_unique<AverageHandler>(*this);
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
    bool as3d = dof.GetNum() == 2; // allows _warp by vector for 2d displacements
    for (size_t iCorner = 0; iCorner < mGeometry.cornerCoordinates.size(); ++iCorner)
    {
        auto dofValues = cell.Interpolate(mGeometry.cornerCoordinates[iCorner], dof);

        if (as3d)
            dofValues = EigenCompanion::To3D(dofValues);

        grid->SetPointData(pointIds[iCorner], dof.GetName(), dofValues);
    }
}

void AverageHandler::CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name, UnstructuredGrid* grid)
{
    grid->DefineCellData(name);
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(values[0].rows(), values[0].cols());
    for (auto val : values)
        sum += val;

    grid->SetCellData(cellId, name, sum / values.size());
}

void AverageHandler::PointData(const CellInterface& cell, std::function<Eigen::VectorXd(Eigen::VectorXd)> f,
                               std::vector<int> pointIds, std::string name, UnstructuredGrid* grid)
{
    grid->DefinePointData(name);
    for (size_t iCorner = 0; iCorner < mGeometry.cornerCoordinates.size(); ++iCorner)
    {
        auto coords = cell.Interpolate(mGeometry.cornerCoordinates[iCorner]);
        auto value = f(coords);
        grid->SetPointData(pointIds[iCorner], name, value);
    }
}
