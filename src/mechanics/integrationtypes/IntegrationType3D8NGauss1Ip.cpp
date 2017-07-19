#include "mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include <cassert>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

//! @brief constructor
NuTo::IntegrationType3D8NGauss1Ip::IntegrationType3D8NGauss1Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType3D8NGauss1Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum == 0);
    return Eigen::Vector3d::Zero();
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType3D8NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType3D8NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 8;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType3D8NGauss1Ip::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                              std::vector<double>& VisualizationPointLocalCoordinates,
                                                              unsigned int& NumVisualizationCells,
                                                              std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                              std::vector<unsigned int>& VisualizationCellsIncidence,
                                                              std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 8;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIP.push_back(0);
}
#endif // ENABLE_VISUALIZE
