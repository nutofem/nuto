#include "visualize/CellGeometryVoronoi.h"
#include "mechanics/interpolation/TypeDefs.h"

namespace NuTo
{
namespace Test
{
template <int TDim>
struct TensorProductGeometryVoronoi
{
    static void buildCellGeometryVoronoi(const std::vector<double>& cellCorners1D,
                                         NuTo::Visualize::CellGeometryVoronoi& cellGeometry)
    {
        std::vector<double> VisualizationPoints1D;
        size_t NumVisualizationPoints1D = cellCorners1D.size() + 1;

        VisualizationPoints1D.push_back(-1.);
        for (size_t i = 1; i < NumVisualizationPoints1D - 1; i++)
        {
            VisualizationPoints1D.push_back(0.5 * (cellCorners1D[i - 1] + cellCorners1D[i]));
        }
        VisualizationPoints1D.push_back(1.);

        size_t NumVisualizationPoints = 1;
        for (size_t dim = 0; dim < TDim; dim++)
            NumVisualizationPoints *= NumVisualizationPoints1D;

        cellGeometry.mCellCornerCoords.reserve(NumVisualizationPoints);

        for (size_t i = 0; i < NumVisualizationPoints; i++)
        {
            NaturalCoords coordinate(TDim);
            size_t power = 1;
            for (size_t dim = 0; dim < TDim; dim++)
            {
                size_t index = (i / power) % NumVisualizationPoints1D;
                coordinate(dim) = VisualizationPoints1D[index];
                power *= NumVisualizationPoints1D;
            }
            cellGeometry.mCellCornerCoords.push_back(coordinate);
        }

        SetCellType(cellGeometry);

        SetVisualizationCells(cellCorners1D.size(), cellGeometry);
    }

    static void SetCellType(NuTo::Visualize::CellGeometryVoronoi& cellGeometry);

    static void SetVisualizationCells(size_t numCellCorners1D, NuTo::Visualize::CellGeometryVoronoi& cellGeometry);
};

template <>
void TensorProductGeometryVoronoi<1>::SetCellType(NuTo::Visualize::CellGeometryVoronoi& cellGeometry)
{
    cellGeometry.mCellType = NuTo::eCellTypes::LINE;
}

template <>
void TensorProductGeometryVoronoi<2>::SetCellType(NuTo::Visualize::CellGeometryVoronoi& cellGeometry)
{
    cellGeometry.mCellType = NuTo::eCellTypes::QUAD;
}

template <>
void TensorProductGeometryVoronoi<3>::SetCellType(NuTo::Visualize::CellGeometryVoronoi& cellGeometry)
{
    cellGeometry.mCellType = NuTo::eCellTypes::HEXAHEDRON;
}

template <>
void TensorProductGeometryVoronoi<1>::SetVisualizationCells(size_t numCellCorners1D,
                                                            NuTo::Visualize::CellGeometryVoronoi& cellGeometry)
{
    for (size_t i = 0; i < numCellCorners1D - 1; i++)
    {
        cellGeometry.mVoronoiCells.push_back({i, i + 1});
    }
}

template <>
void TensorProductGeometryVoronoi<2>::SetVisualizationCells(size_t numCellCorners1D,
                                                            NuTo::Visualize::CellGeometryVoronoi& cellGeometry)
{
    for (size_t row = 0; row < numCellCorners1D; row++)
    {
        for (size_t col = 0; col < numCellCorners1D; col++)
        {
            size_t start = row * (numCellCorners1D + 1) + col;
            cellGeometry.mVoronoiCells.push_back(
                    {start, start + 1, start + numCellCorners1D + 2, start + numCellCorners1D + 1});
        }
    }
}

template <>
void TensorProductGeometryVoronoi<3>::SetVisualizationCells(size_t numCellCorners1D,
                                                            NuTo::Visualize::CellGeometryVoronoi& cellGeometry)
{
    for (size_t height = 0; height < numCellCorners1D; height++)
    {
        for (size_t row = 0; row < numCellCorners1D; row++)
        {
            for (size_t col = 0; col < numCellCorners1D; col++)
            {
                size_t start1 =
                        row * (numCellCorners1D + 1) + col + height * ((numCellCorners1D + 1) * (numCellCorners1D + 1));
                size_t start2 = row * (numCellCorners1D + 1) + col +
                                (height + 1) * ((numCellCorners1D + 1) * (numCellCorners1D + 1));
                cellGeometry.mVoronoiCells.push_back({start1, start1 + 1, start1 + numCellCorners1D + 2,
                                                      start1 + numCellCorners1D + 1, start2, start2 + 1,
                                                      start2 + numCellCorners1D + 2, start2 + numCellCorners1D + 1});
            }
        }
    }
}
}
}
