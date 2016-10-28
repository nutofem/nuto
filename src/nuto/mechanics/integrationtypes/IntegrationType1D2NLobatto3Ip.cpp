#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

NuTo::IntegrationType1D2NLobatto3Ip::IntegrationType1D2NLobatto3Ip():
    iPts{-1.,0.,1.}, weights{0.333333333333333, 1.333333333333333, 0.333333333333333}
{
}


void NuTo::IntegrationType1D2NLobatto3Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates) const
{
    if(rIpNum >= 0 && rIpNum < 3)
        rCoordinates = iPts[rIpNum];
    else
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
}


unsigned int NuTo::IntegrationType1D2NLobatto3Ip::GetNumIntegrationPoints() const
{
    return 3;
}


double NuTo::IntegrationType1D2NLobatto3Ip::GetIntegrationPointWeight(int rIpNum) const
{
    if(rIpNum >= 0 && rIpNum < 3) return weights[rIpNum];
    throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
}


std::string NuTo::IntegrationType1D2NLobatto3Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


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
