#include "TensorProductVoronoiHandler.h"

using namespace NuTo::Visualize;

template <int TDim>
TensorProductVoronoiHandler<TDim>::TensorProductVoronoiHandler(int numCellsPerDirection)
{
    std::vector<double> pointCoordinates1D;
    for (int id = 0; id <= numCellsPerDirection; id++)
        pointCoordinates1D.push_back(-1.0 + 2.0 * id / numCellsPerDirection);

    for (int i = 0; i < std::pow(numCellsPerDirection + 1, TDim); i++)
    {
        Eigen::VectorXd coordinate(TDim);
        int power = 1;
        for (int dim = 0; dim < TDim; dim++)
        {
            int index = (i / power) % (numCellsPerDirection + 1);
            coordinate(dim) = pointCoordinates1D[index];
            power *= (numCellsPerDirection + 1);
        }
        mPointCoordinates.push_back(coordinate);
    }

    SetCellType();

    SetVoronoiCells(numCellsPerDirection);
}

namespace NuTo
{
namespace Visualize
{

template <>
void TensorProductVoronoiHandler<1>::SetCellType()
{
    mCellType = NuTo::eCellTypes::LINE;
}

template <>
void TensorProductVoronoiHandler<2>::SetCellType()
{
    mCellType = NuTo::eCellTypes::QUAD;
}

template <>
void TensorProductVoronoiHandler<3>::SetCellType()
{
    mCellType = NuTo::eCellTypes::HEXAHEDRON;
}

template <>
void TensorProductVoronoiHandler<1>::SetVoronoiCells(int numCellsPerDirection)
{
    for (int i = 0; i < numCellsPerDirection; i++)
    {
        mVoronoiCellCorners.push_back({i, i + 1});
    }
}

template <>
void TensorProductVoronoiHandler<2>::SetVoronoiCells(int numCellsPerDirection)
{
    for (int row = 0; row < numCellsPerDirection; row++)
    {
        for (int col = 0; col < numCellsPerDirection; col++)
        {
            int start = row * (numCellsPerDirection + 1) + col;
            mVoronoiCellCorners.push_back(
                    {start, start + 1, start + numCellsPerDirection + 2, start + numCellsPerDirection + 1});
        }
    }
}

template <>
void TensorProductVoronoiHandler<3>::SetVoronoiCells(int numCellsPerDirection)
{
    for (int height = 0; height < numCellsPerDirection; height++)
    {
        for (int row = 0; row < numCellsPerDirection; row++)
        {
            for (int col = 0; col < numCellsPerDirection; col++)
            {
                int start1 = row * (numCellsPerDirection + 1) + col +
                             height * ((numCellsPerDirection + 1) * (numCellsPerDirection + 1));
                int start2 = row * (numCellsPerDirection + 1) + col +
                             (height + 1) * ((numCellsPerDirection + 1) * (numCellsPerDirection + 1));
                mVoronoiCellCorners.push_back({start1, start1 + 1, start1 + numCellsPerDirection + 2,
                                               start1 + numCellsPerDirection + 1, start2, start2 + 1,
                                               start2 + numCellsPerDirection + 2, start2 + numCellsPerDirection + 1});
            }
        }
    }
}

template <int TDim>
std::vector<int> TensorProductVoronoiHandler<TDim>::WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid)
{
    std::vector<int> pointIds;
    for (auto pointCoords : mPointCoordinates)
        pointIds.push_back(grid->AddPoint(cell.Interpolate(pointCoords)));

    std::vector<int> localSubCells;
    for (auto voronoiCorners : mVoronoiCellCorners)
    {
        std::vector<int> localIds;
        for (int i : voronoiCorners)
            localIds.push_back(pointIds[i]);
        localSubCells.push_back(grid->AddCell(localIds, mCellType));
    }
    mSubCells.push_back(localSubCells);
    return pointIds;
}

template <int TDim>
void TensorProductVoronoiHandler<TDim>::WriteDofValues(const CellInterface& cell, const DofType dof,
                                                       std::vector<int> pointIds, UnstructuredGrid* grid)
{
    grid->DefinePointData(dof.GetName());
    int i = 0;
    for (auto pointCoords : mPointCoordinates)
    {
        grid->SetPointData(pointIds[i], dof.GetName(), cell.Interpolate(pointCoords, dof));
        ++i;
    }
}

template <int TDim>
void TensorProductVoronoiHandler<TDim>::CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name,
                                                 UnstructuredGrid* grid)
{
    grid->DefineCellData(name);
    int i = 0;
    for (auto subCell : mSubCells[cellId])
    {
        grid->SetCellData(subCell, name, values[i]);
        ++i;
    }
}

template class TensorProductVoronoiHandler<1>;
template class TensorProductVoronoiHandler<2>;
template class TensorProductVoronoiHandler<3>;
}
}
