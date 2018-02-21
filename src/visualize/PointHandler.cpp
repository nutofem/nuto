#include "visualize/PointHandler.h"
#include "visualize/UnstructuredGrid.h"
#include "mechanics/cell/CellInterface.h"
#include "math/EigenCompanion.h"

#include <cassert>

using namespace NuTo::Visualize;

PointHandler::PointHandler(std::vector<Eigen::VectorXd> ipCoords)
    : mIntegrationPointCoords(ipCoords)
{
}

PointHandler::PointHandler(const IntegrationTypeBase& integrationType)
{
    for (int i = 0; i < integrationType.GetNumIntegrationPoints(); ++i)
        mIntegrationPointCoords.push_back(integrationType.GetLocalIntegrationPointCoordinates(i));
}

std::unique_ptr<HandlerInterface> PointHandler::Clone() const
{
    return std::make_unique<PointHandler>(*this);
}

std::vector<int> PointHandler::WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid)
{
    std::vector<int> pointIds;
    for (auto ipCoord : mIntegrationPointCoords)
    {
        auto coords = cell.Interpolate(ipCoord);
        pointIds.push_back(grid->AddPoint(coords));
    }
    return pointIds;
}

void PointHandler::WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                  UnstructuredGrid* grid)
{
    grid->DefinePointData(dof.GetName());
    for (size_t i = 0; i < pointIds.size(); ++i)
    {
        auto dofValue = cell.Interpolate(mIntegrationPointCoords[i], dof);
        grid->SetPointData(pointIds[i], dof.GetName(), dofValue);
    }
}

void PointHandler::CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name, UnstructuredGrid* grid)
{
    assert(values.size() == mIntegrationPointCoords.size());
    grid->DefinePointData(name);
    const int offset = cellId * values.size();
    for (size_t i = 0; i < values.size(); ++i)
        grid->SetPointData(i + offset, name, values[i]);
}

void PointHandler::PointData(const CellInterface& cell, std::function<Eigen::VectorXd(Eigen::VectorXd)> f,
                             std::vector<int> pointIds, std::string name, UnstructuredGrid* grid)
{
    grid->DefinePointData(name);
    for (size_t i = 0; i < pointIds.size(); ++i)
    {
        auto coords = cell.Interpolate(mIntegrationPointCoords[i]);
        auto value = f(coords);
        grid->SetPointData(pointIds[i], name, value);
    }
}
