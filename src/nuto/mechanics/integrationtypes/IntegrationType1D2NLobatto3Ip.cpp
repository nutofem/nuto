// $Id: IntegrationType1D2NLobatto3Ip.cpp 345 2010-10-19 07:50:21Z arnold2 $

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include <assert.h>


// constructor
NuTo::IntegrationType1D2NLobatto3Ip::IntegrationType1D2NLobatto3Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType1D2NLobatto3Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    switch (rIpNum)
    {
    case 0 :
        rCoordinates = -1; //
        break;
    case 1 :
        rCoordinates =  0.0;
        break;
    case 2 :
        rCoordinates =  1; //
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NLobatto3Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
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
    switch (rIpNum)
    {
    case 0 :
        return 0.333333333333333; // 1/3
    case 1 :
    	return 1.333333333333333; // 4/3
    case 2 :
        return 0.333333333333333; // 1/3
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NLobatto3Ip::GetIntegrationPointWeight] Ip number out of range.");
    }
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
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 4;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 3;
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
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
