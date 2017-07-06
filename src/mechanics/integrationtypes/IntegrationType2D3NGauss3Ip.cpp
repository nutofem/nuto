#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include <cassert>


//! @brief constructor
NuTo::IntegrationType2D3NGauss3Ip::IntegrationType2D3NGauss3Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss3Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 3);
    switch (rIpNum)
    {
    case 0:
        return Eigen::Vector2d({1. / 6., 1. / 6.});
    case 1:
        return Eigen::Vector2d({4. / 6., 1. / 6.});
    case 2:
        return Eigen::Vector2d({1. / 6., 4. / 6.});
    default:
        throw MechanicsException(
                "[NuTo::IntegrationType2D3NGauss3Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


Eigen::MatrixXd NuTo::IntegrationType2D3NGauss3Ip::GetNaturalIntegrationPointCoordinates() const
{
    Eigen::MatrixXd naturalCoordinates(2, 3);
    naturalCoordinates(0, 0) = 1. / 6.;
    naturalCoordinates(1, 0) = 1. / 6.;

    naturalCoordinates(0, 1) = 4. / 6.;
    naturalCoordinates(1, 1) = 1. / 6.;

    naturalCoordinates(0, 2) = 1. / 6.;
    naturalCoordinates(1, 2) = 4. / 6.;

    return naturalCoordinates;
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss3Ip::GetNumIntegrationPoints() const
{
    return 3;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss3Ip::GetIntegrationPointWeight(int) const
{
    return 1 / 6.;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss3Ip::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                              std::vector<double>& VisualizationPointLocalCoordinates,
                                                              unsigned int& NumVisualizationCells,
                                                              std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                              std::vector<unsigned int>& VisualizationCellsIncidence,
                                                              std::vector<unsigned int>& VisualizationCellsIP) const
{
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
    VisualizationCellsIP.push_back(0);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIP.push_back(1);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(2);
}
#endif // ENABLE_VISUALIZE
