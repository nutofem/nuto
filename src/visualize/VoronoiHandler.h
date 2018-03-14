#pragma once

#include "visualize/HandlerInterface.h"
#include "visualize/VisualizeEnum.h"

namespace NuTo
{
namespace Visualize
{

struct VoronoiCell
{
    std::vector<int> cellCornerIds;
    eCellTypes cellType;
};

struct VoronoiGeometry
{
    std::vector<Eigen::VectorXd> pointCoordinates;
    std::vector<VoronoiCell> voronoiCells;
};

//! Cell handler that subdivides a cell into subcells.
class VoronoiHandler : public HandlerInterface
{
public:
    //! Constructor.
    //! @param geometry definition of arbitrary voronoi cells
    VoronoiHandler(VoronoiGeometry geometry);

    virtual std::unique_ptr<HandlerInterface> Clone() const override;

    virtual std::vector<int> WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid) override;

    virtual void WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                UnstructuredGrid* grid) override;


    virtual void CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name,
                          UnstructuredGrid* grid) override;

    virtual void PointData(const CellInterface& cell, std::function<Eigen::VectorXd(Eigen::VectorXd)> f,
                           std::vector<int> pointIds, std::string name, UnstructuredGrid* grid) override;

private:
    VoronoiGeometry mGeometry;
    std::vector<std::vector<int>> mSubCells;
};

} // namespace Visualize
} // namespace NuTo
