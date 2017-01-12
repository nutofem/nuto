// $Id: IntegrationType1D2NLobatto3Ip.cpp 345 2010-10-19 07:50:21Z arnold2 $

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

// constructor
NuTo::IntegrationType1D2NLobatto3Ip::IntegrationType1D2NLobatto3Ip():
    iPts{-1.,0.,1.}, weights{0.333333333333333, 1.333333333333333, 0.333333333333333}
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType1D2NLobatto3Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    if(rIpNum >= 0 && rIpNum < 3)
        rCoordinates = iPts[rIpNum];
    else
        throw MechanicsException("[NuTo::IntegrationType1D2NLobatto3Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NLobatto3Ip::GetNumIntegrationPoints()const
{
    return 3;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NLobatto3Ip::GetIntegrationPointWeight(int rIpNum)const
{
    if(rIpNum >= 0 && rIpNum < 3) return weights[rIpNum];
    throw MechanicsException("[NuTo::IntegrationType1D2NLobatto3Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NLobatto3Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NLobatto3Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NLOBATTO3IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NLobatto3Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 4;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 3;
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NLobatto3Ip)
#endif
