#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif //ENABLE_VISUALIZE

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss4Ip.h"
#include <assert.h>


NuTo::IntegrationType1D2NGauss4Ip::IntegrationType1D2NGauss4Ip() {}


void NuTo::IntegrationType1D2NGauss4Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates) const
{
    switch (rIpNum)
    {
    case 0 :
        rCoordinates = -0.861136311594052575224;
        break;
    case 1 :
        rCoordinates = -0.339981043584856264803;
        break;
    case 2 :
        rCoordinates = 0.339981043584856264803;
        break;
    case 3 :
        rCoordinates = 0.861136311594052575224;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


unsigned int NuTo::IntegrationType1D2NGauss4Ip::GetNumIntegrationPoints() const
{
    return 4;
}


double NuTo::IntegrationType1D2NGauss4Ip::GetIntegrationPointWeight(int rIpNum) const
{
    switch (rIpNum)
    {
    case 0 :
        return 0.34785484513745385737;
    case 1 :
    	return 0.65214515486254614262;
    case 2 :
    	return 0.65214515486254614262;
    case 3 :
        return 0.34785484513745385737;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


std::string NuTo::IntegrationType1D2NGauss4Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType1D2NGauss4Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NGAUSS4IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NGauss4Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 4;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.3043);
    VisualizationPointLocalCoordinates.push_back(0.3043);
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
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NGauss4Ip)
#endif
