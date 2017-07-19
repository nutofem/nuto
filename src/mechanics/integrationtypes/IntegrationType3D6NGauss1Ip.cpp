#include "mechanics/integrationtypes/IntegrationType3D6NGauss1Ip.h"
#include <cassert>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE


Eigen::VectorXd NuTo::IntegrationType3D6NGauss1Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum == 0);
    return Eigen::Vector3d({1. / 3., 1. / 3., 0});
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D6NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D6NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 1;
}

#ifdef ENABLE_VISUALIZE
NuTo::IntegrationTypeBase::IpCellInfo NuTo::IntegrationType3D6NGauss1Ip::GetVisualizationCells() const
{
    IpCellInfo ipCellInfo;
    constexpr int numPoints = 6;

    ipCellInfo.vertices.resize(numPoints);
    ipCellInfo.vertices[0].localCoords = Eigen::Vector3d(0, 0, -1);
    ipCellInfo.vertices[1].localCoords = Eigen::Vector3d(1, 0, -1);
    ipCellInfo.vertices[2].localCoords = Eigen::Vector3d(0, 1, -1);

    ipCellInfo.vertices[3].localCoords = Eigen::Vector3d(0, 0, 1);
    ipCellInfo.vertices[4].localCoords = Eigen::Vector3d(1, 0, 1);
    ipCellInfo.vertices[5].localCoords = Eigen::Vector3d(0, 1, 1);

    constexpr int numCells = 1;

    ipCellInfo.cells.resize(numCells);
    ipCellInfo.cells[0].cellType = NuTo::eCellTypes::WEDGE;
    ipCellInfo.cells[0].ipId = 0;
    ipCellInfo.cells[0].pointIds = {0, 1, 2, 3, 4, 5};

    return ipCellInfo;
}
#endif // ENABLE_VISUALIZE
