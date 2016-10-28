#include "nuto/mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

NuTo::IntegrationType3D4NGauss1Ip::IntegrationType3D4NGauss1Ip()
{
}


void NuTo::IntegrationType3D4NGauss1Ip::GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3]) const
{
    assert(rIpNum==0);
    rCoordinates[0] = 0.25;
    rCoordinates[1] = 0.25;
    rCoordinates[2] = 0.25;
}


unsigned int NuTo::IntegrationType3D4NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}


double NuTo::IntegrationType3D4NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 1.0/6.0;
}


std::string NuTo::IntegrationType3D4NGauss1Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType3D4NGauss1Ip::GetStrIdentifierStatic()
{
    return std::string("3D4NGAUSS1IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType3D4NGauss1Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 4;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::TETRAEDER);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(0);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType3D4NGauss1Ip)
#endif
