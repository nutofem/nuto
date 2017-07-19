#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include <cassert>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

NuTo::IntegrationType3D4NGauss4Ip::IntegrationType3D4NGauss4Ip()
{
    mCoordinates[0] = Eigen::Vector3d({0.13819660, 0.13819660, 0.13819660});
    mCoordinates[1] = Eigen::Vector3d({0.58541020, 0.13819660, 0.13819660});
    mCoordinates[2] = Eigen::Vector3d({0.13819660, 0.58541020, 0.13819660});
    mCoordinates[3] = Eigen::Vector3d({0.13819660, 0.13819660, 0.58541020});
}

Eigen::VectorXd NuTo::IntegrationType3D4NGauss4Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 4);
    return mCoordinates[rIpNum];
}

int NuTo::IntegrationType3D4NGauss4Ip::GetNumIntegrationPoints() const
{
    return 4;
}

double NuTo::IntegrationType3D4NGauss4Ip::GetIntegrationPointWeight(int) const
{
    return 1. / 24.;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType3D4NGauss4Ip::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                              std::vector<double>& VisualizationPointLocalCoordinates,
                                                              unsigned int& NumVisualizationCells,
                                                              std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                              std::vector<unsigned int>& VisualizationCellsIncidence,
                                                              std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 15; // Identical with the nodes + centroid of faces + centroid

    // Point 0
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 8
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 9
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 10 (centroid of  0 1 3)
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(0.3333333);

    // Point 11 (centroid of 1 2 3)
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.3333333);

    // Point 12 (centoid 0 2 3)
    VisualizationPointLocalCoordinates.push_back(0.);
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.3333333);

    // Point 13 (centroid 0 1 2)
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.3333333);
    VisualizationPointLocalCoordinates.push_back(0.0);

    // Point 14 (centroid)
    VisualizationPointLocalCoordinates.push_back(0.25);
    VisualizationPointLocalCoordinates.push_back(0.25);
    VisualizationPointLocalCoordinates.push_back(0.25);

    NumVisualizationCells = 4;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIP.push_back(0);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIP.push_back(1);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIP.push_back(2);

    // cell 3
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIP.push_back(3);
}
#endif // ENABLE_VISUALIZE
