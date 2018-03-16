#pragma once

#include "nuto/visualize/HandlerInterface.h"
#include "nuto/visualize/VisualizeEnum.h"

namespace NuTo
{
namespace Visualize
{

//! Cell handler that averages the cell values over each cell.
class AverageHandler : public HandlerInterface
{
public:
    virtual std::unique_ptr<HandlerInterface> Clone() const override;

    virtual std::vector<int> WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid) override;

    virtual void WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                UnstructuredGrid* grid) override;

    virtual void CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name,
                          UnstructuredGrid* grid) override;

    virtual void PointData(const CellInterface& cell, std::function<Eigen::VectorXd(Eigen::VectorXd)> f,
                           std::vector<int> pointIds, std::string name, UnstructuredGrid* grid) override;
};

} // namespace Visualize
} // namespace NuTo
