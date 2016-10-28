#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss3Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

NuTo::IntegrationType1D2NGauss3Ip::IntegrationType1D2NGauss3Ip() {}


void NuTo::IntegrationType1D2NGauss3Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates) const
{
    switch (rIpNum)
    {
    case 0 :
        rCoordinates = -0.774596669241483; // -sqr(3/5)
        break;
    case 1 :
        rCoordinates =  0.0;
        break;
    case 2 :
        rCoordinates =  0.774596669241483; // sqr(3/5)
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


unsigned int NuTo::IntegrationType1D2NGauss3Ip::GetNumIntegrationPoints() const
{
    return 3;
}


double NuTo::IntegrationType1D2NGauss3Ip::GetIntegrationPointWeight(int rIpNum) const
{
    switch (rIpNum)
    {
    case 0 :
        return 0.555555555555556; // 5/9
    case 1 :
    	return 0.888888888888889; // 8/9
    case 2 :
        return 0.555555555555556; // 5/9
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


std::string NuTo::IntegrationType1D2NGauss3Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType1D2NGauss3Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NGAUSS3IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NGauss3Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 4;
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.3873);
    VisualizationPointLocalCoordinates.push_back(0.3873);
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
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NGauss3Ip)
#endif
