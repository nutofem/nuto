#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NLobatto.h"


//! @brief constructor
NuTo::IntegrationType2D4NLobatto9Ip::IntegrationType2D4NLobatto9Ip()
{
    NuTo::IntegrationType1D2NLobatto Lobatto1D2N3Ip(3);
    Eigen::Vector3d coordinates1D2N3Ip;
    Eigen::Vector3d weights1D2N3Ip;

    // get the 1D integration point coordinates and weights
    for (int i = 0; i < 3; i++)
    {
        coordinates1D2N3Ip[i] = Lobatto1D2N3Ip.GetLocalIntegrationPointCoordinates(i)[0];
        weights1D2N3Ip[i] = Lobatto1D2N3Ip.GetIntegrationPointWeight(i);
    }

    // calculate the 2D integratration point coordinates and weights
    mWeights.resize(GetNumIntegrationPoints());
    mPts.resize(GetNumIntegrationPoints());
    int ipNum = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            mWeights[ipNum] = weights1D2N3Ip[i] * weights1D2N3Ip[j];
            mPts[ipNum][0] = coordinates1D2N3Ip[j];
            mPts[ipNum][1] = coordinates1D2N3Ip[i];
            ipNum++;
        }
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D4NLobatto9Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if (rIpNum < 0 || rIpNum > 8)
        throw Exception(
                "[NuTo::IntegrationType2D4NLobatto9Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    return mPts[rIpNum];
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D4NLobatto9Ip::GetNumIntegrationPoints() const
{
    return 9;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D4NLobatto9Ip::GetIntegrationPointWeight(int rIpNum) const
{
    if (rIpNum < 0 || rIpNum > 8)
        throw Exception(
                "[NuTo::IntegrationType2D4NLobatto9Ip::GetIntegrationPointWeight] Ip number out of range.");
    return mWeights[rIpNum];
}


#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D4NLobatto9Ip::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                                                std::vector<double>& VisualizationPointLocalCoordinates,
                                                                unsigned int& NumVisualizationCells,
                                                                std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                                                std::vector<unsigned int>& VisualizationCellsIncidence,
                                                                std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 16;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(-0.5);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(1.0);
    VisualizationPointLocalCoordinates.push_back(-1.);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.5);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(-0.5);
    VisualizationPointLocalCoordinates.push_back(-0.5);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(-0.5);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(1.0);
    VisualizationPointLocalCoordinates.push_back(-0.5);

    // Point 8
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 9
    VisualizationPointLocalCoordinates.push_back(-0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 10
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 11
    VisualizationPointLocalCoordinates.push_back(1.0);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 12
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 13
    VisualizationPointLocalCoordinates.push_back(-0.5);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 14
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 15
    VisualizationPointLocalCoordinates.push_back(1.0);
    VisualizationPointLocalCoordinates.push_back(1.);

    NumVisualizationCells = 9;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIP.push_back(0);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIP.push_back(1);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIP.push_back(2);

    // cell 3
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIP.push_back(3);

    // cell 4
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIP.push_back(4);

    // cell 5
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIP.push_back(5);

    // cell 6
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIP.push_back(6);

    // cell 7
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIP.push_back(7);

    // cell 8
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(15);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIP.push_back(8);
}
#endif // ENABLE_VISUALIZE
