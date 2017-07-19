#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include <cassert>


//! @brief constructor
NuTo::IntegrationType2D4NGauss1Ip::IntegrationType2D4NGauss1Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D4NGauss1Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum == 0);
    return Eigen::Vector2d({0., 0.});
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D4NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D4NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 4;
}


#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D4NGauss1Ip::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                              std::vector<double>& VisualizationPointLocalCoordinates,
                                                              unsigned int& NumVisualizationCells,
                                                              std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                              std::vector<unsigned int>& VisualizationCellsIncidence,
                                                              std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 4;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(0);
}
#endif // ENABLE_VISUALIZE
