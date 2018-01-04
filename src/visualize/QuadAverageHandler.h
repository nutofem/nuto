#pragma once

#include "mechanics/cell/CellInterface.h"
#include "visualize/UnstructuredGrid.h"

namespace NuTo
{
namespace Visualize
{

class QuadAverageHandler
{
public:
    std::vector<int> WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid);

    void WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                        UnstructuredGrid* grid);

    void CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name, UnstructuredGrid* grid);

private:
    const std::vector<Eigen::VectorXd> cornerCoordinates = {
            Eigen::Vector3d(-1.0, -1.0, 0.0), Eigen::Vector3d(1.0, -1.0, 0.0), Eigen::Vector3d(1.0, 1.0, 0.0),
            Eigen::Vector3d(-1.0, 1.0, 0.0)};
};

} // namespace Visualize
} // namespace NuTo
