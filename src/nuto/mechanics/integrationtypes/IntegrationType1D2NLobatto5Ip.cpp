#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

NuTo::IntegrationType1D2NLobatto5Ip::IntegrationType1D2NLobatto5Ip():
    iPts{-1., -0.654653670707977087, 0., +0.654653670707977087, +1.},
    weights{0.1, 0.544444444444444, 0.711111111111111, 0.544444444444444, 0.1}
{
}

void NuTo::IntegrationType1D2NLobatto5Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates) const
{
    if(rIpNum >= 0 && rIpNum < 5)
        rCoordinates = iPts[rIpNum];
    else
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
}


unsigned int NuTo::IntegrationType1D2NLobatto5Ip::GetNumIntegrationPoints() const
{
    return 5;
}

double NuTo::IntegrationType1D2NLobatto5Ip::GetIntegrationPointWeight(int rIpNum) const
{
    if(rIpNum >= 0 && rIpNum < 5) return weights[rIpNum];
    throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
}


std::string NuTo::IntegrationType1D2NLobatto5Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType1D2NLobatto5Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NLOBATTO5IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NLobatto5Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
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
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
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
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
    VisualizationCellsIP.push_back(3);
    VisualizationCellsIP.push_back(4);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NLobatto5Ip)
#endif
