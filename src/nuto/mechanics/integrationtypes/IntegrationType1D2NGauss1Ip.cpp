#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include <assert.h>

NuTo::IntegrationType1D2NGauss1Ip::IntegrationType1D2NGauss1Ip() {}

void NuTo::IntegrationType1D2NGauss1Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates) const
{
    assert(rIpNum == 0);
    rCoordinates = 0;
}


unsigned int NuTo::IntegrationType1D2NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}


double NuTo::IntegrationType1D2NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 2;
}


std::string NuTo::IntegrationType1D2NGauss1Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType1D2NGauss1Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NGAUSS1IP");
}


#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NGauss1Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 2;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 1;
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIP.push_back(0);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NGauss1Ip)
#endif
