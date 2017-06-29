#include "mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

// constructor
NuTo::IntegrationType1D2NLobatto4Ip::IntegrationType1D2NLobatto4Ip():
    iPts{-1., -0.447213595499957928, 0.447213595499957928, 1.},
    weights{0.16666666666666667, 0.83333333333333333, 0.83333333333333333, 0.16666666666666667}
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType1D2NLobatto4Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if(rIpNum >= 0 && rIpNum < 4)
        return Eigen::Matrix<double, 1, 1>::Constant(iPts[rIpNum]);
    else
        throw Exception("[NuTo::IntegrationType1D2NLobatto4Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NLobatto4Ip::GetNumIntegrationPoints()const
{
    return 4;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NLobatto4Ip::GetIntegrationPointWeight(int rIpNum)const
{
    if(rIpNum >= 0 && rIpNum < 4) return weights[rIpNum];
    throw Exception("[NuTo::IntegrationType1D2NLobatto4Ip::GetIntegrationPointWeight] Ip number out of range.");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NLobatto4Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 5;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.7236);
    VisualizationPointLocalCoordinates.push_back(0.);
    VisualizationPointLocalCoordinates.push_back(0.7236);
    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 4;
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
    VisualizationCellsIP.push_back(3);
}
#endif // ENABLE_VISUALIZE


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NLobatto4Ip)
#endif
