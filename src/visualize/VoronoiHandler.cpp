#include "visualize/VoronoiHandler.h"
#include "visualize/UnstructuredGrid.h"
#include "mechanics/cell/CellInterface.h"
#include "math/EigenCompanion.h"
#include "visualize/VoronoiGeometries.h"

using namespace NuTo;
using namespace NuTo::Visualize;

std::unique_ptr<HandlerInterface> VoronoiHandler::Clone() const
{
    return std::make_unique<VoronoiHandler>(*this);
}

int GetNumSubCells(const CellInterface& cell)
{
    auto nullFct = [](const CellIpData&) { return Eigen::Matrix<double, 1, 1>(); };
    auto results = cell.Eval(nullFct);
    return results.size();
}

VoronoiGeometry CreateGeometry(const CellInterface& cell)
{
    switch (cell.GetShape().Enum())
    {
    case eShape::Line:
        return VoronoiGeometryLine(GetNumSubCells(cell));
    case eShape::Quadrilateral:
        return VoronoiGeometryQuad(std::sqrt(GetNumSubCells(cell)));
    case eShape::Hexahedron:
        return VoronoiGeometryBrick(std::cbrt(GetNumSubCells(cell)));
    default:
        throw Exception(__PRETTY_FUNCTION__, "No Voronoi Geometry defined for this shape.");
    }
}

VoronoiGeometry VoronoiHandler::GetGeometry(const CellInterface& cell)
{
    auto shapeNumSubCellsPair = std::make_pair(cell.GetShape().Enum(), GetNumSubCells(cell));
    auto search = mGeometries.find(shapeNumSubCellsPair);
    if (search != mGeometries.end())
        return search->second;
    else
    {
        mGeometries[shapeNumSubCellsPair] = CreateGeometry(cell);
        return mGeometries[shapeNumSubCellsPair];
    }
}

std::vector<int> VoronoiHandler::WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid)
{
    std::vector<int> pointIds;
    VoronoiGeometry geometry = GetGeometry(cell);
    for (auto pointCoords : geometry.pointCoordinates)
        pointIds.push_back(grid->AddPoint(cell.Interpolate(pointCoords)));

    std::vector<int> localSubCells;
    for (auto voronoiCell : geometry.voronoiCells)
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
    bool as3d = dof.GetNum() == 2; // allows _warp by vector for 2d displacements
    VoronoiGeometry geometry = GetGeometry(cell);
    for (size_t iPoint = 0; iPoint < geometry.pointCoordinates.size(); ++iPoint)
    {
        auto dofValues = cell.Interpolate(geometry.pointCoordinates[iPoint], dof);

        if (as3d)
            dofValues = EigenCompanion::To3D(dofValues);

        grid->SetPointData(pointIds[iPoint], dof.GetName(), dofValues);
    }
}

void VoronoiHandler::CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name, UnstructuredGrid* grid)
{
    assert(values.size() == mSubCells[cellId].size());
    grid->DefineCellData(name);
    for (size_t iValue = 0; iValue < values.size(); ++iValue)
        grid->SetCellData(mSubCells[cellId][iValue], name, values[iValue]);
}

void VoronoiHandler::PointData(const CellInterface& cell, std::function<Eigen::VectorXd(Eigen::VectorXd)> f,
                               std::vector<int> pointIds, std::string name, UnstructuredGrid* grid)
{
    grid->DefinePointData(name);
    VoronoiGeometry geometry = GetGeometry(cell);
    for (size_t iPoint = 0; iPoint < geometry.pointCoordinates.size(); ++iPoint)
    {
        auto coords = cell.Interpolate(geometry.pointCoordinates[iPoint]);
        auto value = f(coords);
        grid->SetPointData(pointIds[iPoint], name, value);
    }
}
