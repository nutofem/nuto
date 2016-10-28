#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

NuTo::IntegrationType1D2NGauss2Ip::IntegrationType1D2NGauss2Ip() {}


void NuTo::IntegrationType1D2NGauss2Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates) const
{
    assert(rIpNum==0 || rIpNum==1);
    switch (rIpNum)
    {
    case 0 :
        rCoordinates = -0.577350269189626;
        break;
    case 1 :
        rCoordinates =  0.577350269189626;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


unsigned int NuTo::IntegrationType1D2NGauss2Ip::GetNumIntegrationPoints() const
{
    return 2;
}


double NuTo::IntegrationType1D2NGauss2Ip::GetIntegrationPointWeight(int) const
{
    return 1.0;
}


std::string NuTo::IntegrationType1D2NGauss2Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType1D2NGauss2Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NGAUSS2IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NGauss2Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 3;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 2;
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NGauss2Ip)
#endif
