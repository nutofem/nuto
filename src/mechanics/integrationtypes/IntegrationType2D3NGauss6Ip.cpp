#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"
#include <cassert>

//! @brief constructor
NuTo::IntegrationType2D3NGauss6Ip::IntegrationType2D3NGauss6Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss6Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 6);

    const double a = 0.445948490915965;
    const double b = 0.091576213509771;

    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector2d({a, a});
    case 1:
        return Eigen::Vector2d({1 - 2 * a, a});
    case 2:
        return Eigen::Vector2d({a, 1 - 2 * a});
    case 3:
        return Eigen::Vector2d({b, b});
    case 4:
        return Eigen::Vector2d({1 - 2 * b, b});
    case 5:
        return Eigen::Vector2d({b, 1 - 2 * b});
    default:
        throw MechanicsException(
                "[NuTo::IntegrationType2D3NGauss6Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss6Ip::GetNumIntegrationPoints() const
{
    return 6;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss6Ip::GetIntegrationPointWeight(int rIpNum) const
{
    const double c = 0.111690794839005;
    const double d = 0.054975871827661;

    assert(rIpNum >= 0 && rIpNum < 6);
    switch (rIpNum)
    {
    case 0:
        return c;
    case 1:
        return c;
    case 2:
        return c;
    case 3:
        return d;
    case 4:
        return d;
    case 5:
        return d;
    default:
        throw MechanicsException(
                "[NuTo::IntegrationType2D3NGauss6Ip::GetIntegrationPointWeight] Ip number out of range.");
    }
}


#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss6Ip::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                              std::vector<double>& VisualizationPointLocalCoordinates,
                                                              unsigned int& NumVisualizationCells,
                                                              std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                              std::vector<unsigned int>& VisualizationCellsIncidence,
                                                              std::vector<unsigned int>& VisualizationCellsIP) const
{

    // only 3 integration points (1,2,3) are visualised. TODO: Voronoi decomposition + triangulation for proper
    // visualisation

    NumVisualizationPoints = 7;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(1. / 3.);
    VisualizationPointLocalCoordinates.push_back(1. / 3.);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    NumVisualizationCells = 3;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(3);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIP.push_back(4);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(5);
}
#endif // ENABLE_VISUALIZE
