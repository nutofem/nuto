#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/UnstructuredGrid.h"
#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/math/EigenCompanion.h"
#include "nuto/visualize/AverageGeometries.h"

using namespace NuTo;
using namespace NuTo::Visualize;

std::unique_ptr<HandlerInterface> AverageHandler::Clone() const
{
    return std::make_unique<AverageHandler>(*this);
}

AverageGeometry GetGeometry(const Shape& shape)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return AverageGeometryLine();
    case eShape::Quadrilateral:
        return AverageGeometryQuad();
    case eShape::Triangle:
        return AverageGeometryTriangle();
    case eShape::Tetrahedron:
        return AverageGeometryTetrahedron();
    case eShape::Prism:
        return AverageGeometryPrism();
    default:
        throw Exception(__PRETTY_FUNCTION__, "No AverageGeometry defined for this shape.");
    }
}
std::vector<int> AverageHandler::WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid)
{
    std::vector<int> pointIds;
    AverageGeometry geometry = GetGeometry(cell.GetShape());
    pointIds.reserve(geometry.cornerCoordinates.size());

    for (auto corner : geometry.cornerCoordinates)
        pointIds.push_back(grid->AddPoint(cell.Interpolate(corner)));
    grid->AddCell(pointIds, geometry.cellType);
    return pointIds;
}

void AverageHandler::WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                    UnstructuredGrid* grid)
{
    grid->DefinePointData(dof.GetName());
    bool as3d = dof.GetNum() == 2; // allows _warp by vector for 2d displacements
    AverageGeometry geometry = GetGeometry(cell.GetShape());
    for (size_t iCorner = 0; iCorner < geometry.cornerCoordinates.size(); ++iCorner)
    {
        auto dofValues = cell.Interpolate(geometry.cornerCoordinates[iCorner], dof);

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
    AverageGeometry geometry = GetGeometry(cell.GetShape());
    for (size_t iCorner = 0; iCorner < geometry.cornerCoordinates.size(); ++iCorner)
    {
        auto coords = cell.Interpolate(geometry.cornerCoordinates[iCorner]);
        auto value = f(coords);
        grid->SetPointData(pointIds[iCorner], name, value);
    }
}
