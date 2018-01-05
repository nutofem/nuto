#pragma once

#include "mechanics/cell/CellInterface.h"
#include "visualize/UnstructuredGrid.h"

namespace NuTo
{
namespace Visualize
{

template <int TDim>
class TensorProductVoronoiHandler
{
public:
    TensorProductVoronoiHandler(int numCellsPerDirection);

    std::vector<int> WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid);

    void WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                        UnstructuredGrid* grid);

    void CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name, UnstructuredGrid* grid);

private:
    void SetCellType();
    void SetVoronoiCells(int numCellsPerDirection);

    std::vector<std::vector<int>> mVoronoiCellCorners;

    //! for each "mechanics cell", save the IDs of the "visualize cells"
    std::vector<std::vector<int>> mSubCells;
    std::vector<Eigen::VectorXd> mPointCoordinates;
    eCellTypes mCellType;
};

} // namespace Visualize
} // namespace NuTo
