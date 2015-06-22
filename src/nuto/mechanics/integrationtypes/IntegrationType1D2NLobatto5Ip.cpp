// $Id: IntegrationType1D2NLobatto5Ip.cpp 345 2010-10-19 07:50:21Z arnold2 $

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include <assert.h>


// constructor
NuTo::IntegrationType1D2NLobatto5Ip::IntegrationType1D2NLobatto5Ip():
    iPts{-1., -0.654653670707977087, 0., +0.654653670707977087, +1.},
    weights{0.1, 0.544444444444444, 0.711111111111111, 0.544444444444444, 0.1}
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType1D2NLobatto5Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    if(rIpNum >= 0 && rIpNum < 5)
        rCoordinates = iPts[rIpNum];
    else
        throw MechanicsException("[NuTo::IntegrationType1D2NLobatto5Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NLobatto5Ip::GetNumIntegrationPoints()const
{
    return 5;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NLobatto5Ip::GetIntegrationPointWeight(int rIpNum)const
{
    if(rIpNum >= 0 && rIpNum < 5) return weights[rIpNum];
    throw MechanicsException("[NuTo::IntegrationType1D2NLobatto5Ip::GetIntegrationPointWeight] Ip number out of range.");
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NLobatto5Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NLobatto5Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NLOBATTO5IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NLobatto5Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 6;
    VisualizationPointLocalCoordinates.push_back(-1.);
    VisualizationPointLocalCoordinates.push_back(-0.8273);
    VisualizationPointLocalCoordinates.push_back(-0.273);
    VisualizationPointLocalCoordinates.push_back( 0.273);
    VisualizationPointLocalCoordinates.push_back( 0.8273);
    VisualizationPointLocalCoordinates.push_back( 1.);
    NumVisualizationCells = 5;
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
    VisualizationCellsIP.push_back(3);
    VisualizationCellsIP.push_back(4);
}
#endif // ENABLE_VISUALIZE
