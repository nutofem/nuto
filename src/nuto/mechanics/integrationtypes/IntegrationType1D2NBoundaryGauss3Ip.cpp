#include "nuto/mechanics/integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

NuTo::IntegrationType1D2NBoundaryGauss3Ip::IntegrationType1D2NBoundaryGauss3Ip() {}


void NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates) const
{
    switch (rIpNum)
    {
    case 0 :
        rCoordinates = -1; // this is the one located on the boundary of the real boundary element
        break;
    case 1 :
        rCoordinates = -0.774596669241483; // -sqr(3/5)
        break;
    case 2 :
        rCoordinates =  0.0;
        break;
    case 3 :
        rCoordinates =  0.774596669241483; // sqr(3/5)
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


unsigned int NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetNumIntegrationPoints() const
{
    return 4;
}


double NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetIntegrationPointWeight(int rIpNum) const
{
    switch (rIpNum)
    {
    case 0 :
        return 0.0; // 5/9
    case 1 :
        return 0.555555555555556; // 5/9
    case 2 :
    	return 0.888888888888889; // 8/9
    case 3 :
        return 0.555555555555556; // 5/9
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


std::string NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetStrIdentifierStatic()
{
    return std::string("1D2NBOUNDARYGAUSS3IP");
}


#ifdef ENABLE_SERIALIZATION
    template<class Archive>
    void NuTo::IntegrationType1D2NBoundaryGauss3Ip::serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize IntegrationType1D2NBoundaryGauss3Ip" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuTo::IntegrationType1D);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize IntegrationType1D2NBoundaryGauss3Ip" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetVisualizationCells(
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
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
    VisualizationCellsIP.push_back(3);
}
#endif // ENABLE_VISUALIZE
