// $Id$

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss5Ip.h"
#include <assert.h>


// constructor
NuTo::IntegrationType1D2NGauss5Ip::IntegrationType1D2NGauss5Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType1D2NGauss5Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    switch (rIpNum)
    {
    case 0 :
        rCoordinates = -0.906179845938664;
        break;
    case 1 :
        rCoordinates = -0.538469310105683;
        break;
    case 2 :
        rCoordinates = 0.0;
        break;
    case 3 :
        rCoordinates = 0.538469310105683;
        break;
    case 4 :
        rCoordinates = 0.906179845938664;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NGauss5Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NGauss5Ip::GetNumIntegrationPoints()const
{
    return 5;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NGauss5Ip::GetIntegrationPointWeight(int rIpNum)const
{
    switch (rIpNum)
    {
    case 0 :
        return 0.236926885056189;
    case 1 :
    	return 0.478628670499366;
    case 2 :
    	return 0.568888888888889;
    case 3 :
        return 0.478628670499366;
    case 4 :
        return 0.236926885056189;
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NGauss5Ip::GetIntegrationPointWeight] Ip number out of range.");
    }
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NGauss5Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NGauss5Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NGAUSS5IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NGauss5Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 5;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.44444);
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(0.444444);
    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 4;
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
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
    VisualizationCellsIP.push_back(3);
}
#endif // ENABLE_VISUALIZE
