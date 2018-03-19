#pragma once

#include "nuto/visualize/HandlerInterface.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
namespace Visualize
{
//! Visualizes cell and point data only at the integration points (no visualization cells)
class PointHandler : public HandlerInterface
{
public:
    //! ctor
    //! @param ipCoords vector of integration point coordinates
    PointHandler(std::vector<Eigen::VectorXd> ipCoords);

    //! ctor
    //! @param integrationType integration type
    PointHandler(const IntegrationTypeBase& integrationType);

    virtual std::unique_ptr<HandlerInterface> Clone() const override;

    virtual std::vector<int> WriteGeometry(const CellInterface& cell, UnstructuredGrid* grid) override;

    virtual void WriteDofValues(const CellInterface& cell, const DofType dof, std::vector<int> pointIds,
                                UnstructuredGrid* grid) override;

    virtual void CellData(int cellId, std::vector<Eigen::VectorXd> values, std::string name,
                          UnstructuredGrid* grid) override;

    virtual void PointData(const CellInterface& cell, std::function<Eigen::VectorXd(Eigen::VectorXd)> f,
                           std::vector<int> pointIds, std::string name, UnstructuredGrid* grid) override;

private:
    std::vector<Eigen::VectorXd> mIntegrationPointCoords;
};

} /* Visualize */
} /* NuTo */
