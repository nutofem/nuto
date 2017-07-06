#include "mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.h"
#include <cassert>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

//! @brief constructor
NuTo::IntegrationType3D6NGauss2x3Ip::IntegrationType3D6NGauss2x3Ip()
{
}

Eigen::VectorXd NuTo::IntegrationType3D6NGauss2x3Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    /**
     *  x and y [index 0 and 1] taken from 2D3N3IP triangle
     *  z taken from 1D2N2IP truss
     */

    assert(rIpNum >= 0 && rIpNum < 6);
    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector3d({1. / 6., 1. / 6., -0.577350269189626});
    case 1:
        return Eigen::Vector3d({4. / 6., 1. / 6., -0.577350269189626});
    case 2:
        return Eigen::Vector3d({1. / 6., 4. / 6., -0.577350269189626});
    case 3:
        return Eigen::Vector3d({1. / 6., 1. / 6., 0.577350269189626});
    case 4:
        return Eigen::Vector3d({4. / 6., 1. / 6., 0.577350269189626});
    case 5:
        return Eigen::Vector3d({1. / 6., 4. / 6., 0.577350269189626});
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D6NGauss2x3Ip::GetNumIntegrationPoints() const
{
    return 6;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D6NGauss2x3Ip::GetIntegrationPointWeight(int) const
{
    return 1 / 6.;
}

#ifdef ENABLE_VISUALIZE
NuTo::IntegrationTypeBase::IpCellInfo NuTo::IntegrationType3D6NGauss2x3Ip::GetVisualizationCells() const
{
    IpCellInfo ipCellInfo;
    constexpr int numPoints = 21;

    ipCellInfo.vertices.resize(numPoints);
    ipCellInfo.vertices[0].localCoords = Eigen::Vector3d(0, 0, -1);
    ipCellInfo.vertices[1].localCoords = Eigen::Vector3d(0.5, 0, -1);
    ipCellInfo.vertices[2].localCoords = Eigen::Vector3d(1, 0, -1);
    ipCellInfo.vertices[3].localCoords = Eigen::Vector3d(0, 0.5, -1);
    ipCellInfo.vertices[4].localCoords = Eigen::Vector3d(1. / 3., 1. / 3., -1);
    ipCellInfo.vertices[5].localCoords = Eigen::Vector3d(0.5, 0.5, -1);
    ipCellInfo.vertices[6].localCoords = Eigen::Vector3d(0, 1, -1);

    ipCellInfo.vertices[7].localCoords = Eigen::Vector3d(0, 0, 0);
    ipCellInfo.vertices[8].localCoords = Eigen::Vector3d(0.5, 0, 0);
    ipCellInfo.vertices[9].localCoords = Eigen::Vector3d(1, 0, 0);
    ipCellInfo.vertices[10].localCoords = Eigen::Vector3d(0, 0.5, 0);
    ipCellInfo.vertices[11].localCoords = Eigen::Vector3d(1. / 3., 1. / 3., 0);
    ipCellInfo.vertices[12].localCoords = Eigen::Vector3d(0.5, 0.5, 0);
    ipCellInfo.vertices[13].localCoords = Eigen::Vector3d(0, 1, 0);

    ipCellInfo.vertices[14].localCoords = Eigen::Vector3d(0, 0, 1);
    ipCellInfo.vertices[15].localCoords = Eigen::Vector3d(0.5, 0, 1);
    ipCellInfo.vertices[16].localCoords = Eigen::Vector3d(1, 0, 1);
    ipCellInfo.vertices[17].localCoords = Eigen::Vector3d(0, 0.5, 1);
    ipCellInfo.vertices[18].localCoords = Eigen::Vector3d(1. / 3., 1. / 3., 1);
    ipCellInfo.vertices[19].localCoords = Eigen::Vector3d(0.5, 0.5, 1);
    ipCellInfo.vertices[20].localCoords = Eigen::Vector3d(0, 1, 1);

    constexpr int numCells = 6;

    ipCellInfo.cells.resize(numCells);
    ipCellInfo.cells[0].cellType = NuTo::eCellTypes::HEXAHEDRON;
    ipCellInfo.cells[0].ipId = 0;
    ipCellInfo.cells[0].pointIds = {0, 1, 4, 3, 7, 8, 11, 10};

    ipCellInfo.cells[1].cellType = NuTo::eCellTypes::HEXAHEDRON;
    ipCellInfo.cells[1].ipId = 1;
    ipCellInfo.cells[1].pointIds = {1, 2, 5, 4, 8, 9, 12, 11};

    ipCellInfo.cells[2].cellType = NuTo::eCellTypes::HEXAHEDRON;
    ipCellInfo.cells[2].ipId = 2;
    ipCellInfo.cells[2].pointIds = {3, 4, 5, 6, 10, 11, 12, 13};

    ipCellInfo.cells[3].cellType = NuTo::eCellTypes::HEXAHEDRON;
    ipCellInfo.cells[3].ipId = 3;
    ipCellInfo.cells[3].pointIds = {7, 8, 11, 10, 14, 15, 18, 17};

    ipCellInfo.cells[4].cellType = NuTo::eCellTypes::HEXAHEDRON;
    ipCellInfo.cells[4].ipId = 4;
    ipCellInfo.cells[4].pointIds = {8, 9, 12, 11, 15, 16, 19, 18};

    ipCellInfo.cells[5].cellType = NuTo::eCellTypes::HEXAHEDRON;
    ipCellInfo.cells[5].ipId = 5;
    ipCellInfo.cells[5].pointIds = {10, 11, 12, 13, 17, 18, 19, 20};

    return ipCellInfo;
}
#endif // ENABLE_VISUALIZE
