#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

NuTo::IntegrationType3D8NGauss1Ip::IntegrationType3D8NGauss1Ip() {}

void NuTo::IntegrationType3D8NGauss1Ip::GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3]) const
{
    assert(rIpNum==0);
	rCoordinates[0] = 0.;
	rCoordinates[1] = 0.;
	rCoordinates[2] = 0.;
}


unsigned int NuTo::IntegrationType3D8NGauss1Ip::GetNumIntegrationPoints() const
{
    return 1;
}


double NuTo::IntegrationType3D8NGauss1Ip::GetIntegrationPointWeight(int) const
{
    return 8.0;
}


std::string NuTo::IntegrationType3D8NGauss1Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType3D8NGauss1Ip::GetStrIdentifierStatic()
{
    return std::string("3D8NGAUSS1IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType3D8NGauss1Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 8;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIP.push_back(0);
}
#endif // ENABLE_VISUALIZE
