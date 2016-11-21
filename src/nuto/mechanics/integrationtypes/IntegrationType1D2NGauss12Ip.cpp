// $Id$

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss12Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

// constructor
NuTo::IntegrationType1D2NGauss12Ip::IntegrationType1D2NGauss12Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType1D2NGauss12Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    switch (rIpNum)
    {
    case 0 :
        rCoordinates = -0.981560634246719;
        break;
    case 1 :
        rCoordinates = -0.904117256370475;
        break;
    case 2 :
        rCoordinates = -0.769902674194305;
        break;
    case 3 :
        rCoordinates = -0.587317954286617;
        break;
    case 4 :
        rCoordinates = -0.367831498998180;
        break;
    case 5 :
        rCoordinates = -0.125233408511469;
        break;
    case 6 :
        rCoordinates = 0.125233408511469;
        break;
    case 7 :
        rCoordinates = 0.367831498998180;
        break;
    case 8 :
        rCoordinates = 0.587317954286617;
        break;
    case 9 :
        rCoordinates = 0.769902674194305;
        break;
    case 10 :
        rCoordinates = 0.904117256370475;
        break;
    case 11 :
        rCoordinates = 0.981560634246719;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NGauss12Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NGauss12Ip::GetNumIntegrationPoints()const
{
    return 12;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NGauss12Ip::GetIntegrationPointWeight(int rIpNum)const
{
    switch (rIpNum)
    {
    case 0 :
        return 0.047175336386512;
    case 1 :
        return 0.106939325995318;
    case 2 :
        return  0.160078328543346;
    case 3 :
        return 0.203167426723066;
    case 4 :
        return 0.233492536538355;
    case 5 :
        return  0.249147045813403;
    case 6 :
        return 0.249147045813403;
    case 7 :
        return 0.233492536538355;
    case 8 :
        return 0.203167426723066;
    case 9 :
        return 0.160078328543346;
    case 10 :
        return 0.106939325995318;
    case 11 :
        return 0.047175336386512;
    default:
        throw MechanicsException("[NuTo::IntegrationType1D2NGauss12Ip::GetIntegrationPointWeight] Ip number out of range.");
    }
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NGauss12Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NGauss12Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NGAUSS12IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NGauss12Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    double coordinates[2];
    NumVisualizationPoints = GetNumIntegrationPoints()+1;
    VisualizationPointLocalCoordinates.push_back(-1.);

    for(int i = 1; i < GetNumIntegrationPoints(); i++)
    {
        GetLocalIntegrationPointCoordinates1D(i-1, coordinates[0]);
        GetLocalIntegrationPointCoordinates1D(i, coordinates[1]);
        VisualizationPointLocalCoordinates.push_back((coordinates[0] + coordinates[1])/2);
    }
    VisualizationPointLocalCoordinates.push_back(1.);

    NumVisualizationCells = GetNumIntegrationPoints();
    for(int i = 0; i < GetNumIntegrationPoints(); i++)  VisualizationCellType.push_back(NuTo::eCellTypes::LINE);

    VisualizationCellsIncidence.push_back(0);
    for(int i = 1; i < GetNumIntegrationPoints(); i++)
        for(int j = 0; j < 2; j++) VisualizationCellsIncidence.push_back(i);
    VisualizationCellsIncidence.push_back(12);

    for(int i = 0; i < GetNumIntegrationPoints(); i++) VisualizationCellsIP.push_back(i);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NGauss12Ip)
#endif
